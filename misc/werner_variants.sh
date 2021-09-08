GATK=/home/werner/tools/gatk4.1.3.0/gatk/gatk
igvtools=/home/werner/tools/IGVTools/igvtools.jar
samtools=/home/werner/tools/samtools-1.9/samtools
genome_fa=/home/werner/xchrom_snp_skew/data/annot_genomes/GRCh38_Gencode25/GRCh38.p7.genome.fa #Genome fasta file location


inDir=$1  #Specifies BAM file location
chr=$2    #Specifies which chromosome to process
outPrefix=`pwd`
echo $outPrefix
inI=$((SGE_TASK_ID-1))
echo $inI
outDir=$inDir
cd $outDir

chr_fa=/home/werner/xchrom_snp_skew/data/annot_genomes/GRCh38_Gencode25/$chr.fa   #Chromosome fasta file location

echo "Preparing BAM file for SNP calling"

if [ ! -e "$inDir/Aligned.sortedByCoord.out.bam" ]       #Sort the bam if not already done so
then
	$samtools sort -@ 10 $inDir/Aligned.out.bam Aligned.sortedByCoord.out
fi

if [ ! -e "$inDir/Aligned.sortedByCoord.out.bai" ]       #Generate the indexed BAM if not already done so
then
	$samtools index $inDir/Aligned.sortedByCoord.out.bam
fi

$samtools view -b $inDir/Aligned.sortedByCoord.out.bam $chr > $outDir/$chr.bam    #Extract chromosome bam
$samtools view -b -q 10 $outDir/$chr.bam > $outDir/$chr.filt.bam                  #Filter chromosome bam?

echo "Adding read groups to bam file"
$GATK AddOrReplaceReadGroups \
    -I=$outDir/$chr.filt.bam \
    -O=$outDir/$chr.rg.bam \
    -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample

echo "Marking duplicates"
$GATK MarkDuplicates \
    -I=$outDir/$chr.rg.bam  \
    -O=$outDir/$chr.dedupped.bam \
    -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics

echo "Splitting and trimming"
$GATK SplitNCigarReads         \
   -R $genome_fa \
   -I $outDir/$chr.dedupped.bam \
   -O $outDir/$chr.split.filtered.bam 

#Generates the initial .vcf file
echo "Haplotype calling (?)"
$GATK HaplotypeCaller \
   -R $genome_fa \
   -L $chr \
   -I $outDir/$chr.split.filtered.bam \
   --dont-use-soft-clipped-bases -stand-call-conf 0.0 \
   -O $outDir/$chr.filtered.vcf

#Generates the .wig file (contains the reads per position we care about)
echo "Counting SNPs"
java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   $outDir/$chr.split.filtered.bam   \
   $outDir/$chr.split.filtered.wig   \
   $chr_fa

#Grab the SNPs from the .vcf file
#Just grab the SNPs, not worrying about the INDELs
$GATK SelectVariants \
        -R $genome_fa \
        -L $chr \
        -V $chr.filtered.vcf \
        -O $chr.snps.vcf \
        -select-type SNP

#Filter on certain thresholds. These are the recommended annotations and thresholds to filter from GATK
#Only filtering on the SOR, FS, and ReadPosRankSum tests downstream, but good to include the others
$GATK VariantFiltration \
        -R $genome_fa \
        -L $chr \
        -V $chr.snps.vcf \
        -O $chr.snps_filtered.vcf \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

#Convert to a tab delimited file for easy load up in R
#These will be all the SNPs in a sample, with each variant's associated statistics
#The chrX.filtered.vcf will contain SNPs and INDELs that pass dbSNP filtering (next step, cut_files.sh script)
#the chrX.cut and the chrX.snps_filtered.tab (output below) will need to be matched up, done in R 
$GATK VariantsToTable \
        -V $chr.snps_filtered.vcf \
        -O $chr.snps_filtered.tab \
        -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F FS -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR -GF GT -GF AD -GF DP -GF GQ -GF PL \
        --show-filtered true


rm $outDir/$chr.bam $outDir/$chr.filt.bam $outDir/$chr.rg.bam $outDir/$chr.dedupped.bam $outDir/$chr.dedupped.bai $outDir/$chr.split.filtered.bam  $outDir/$chr.split.filtered.bai $outDir/output.metrics $outDir/igv.log $outDir/Aligned.sortedByCoord.out.bam.bai
