#Command line examples for processing GTEx data

#Grab some sample names, I've been processing them in batches of 100.
#For this data, typical sample name looks like this GTEX-ZPIC-0526-SM-DO91W
files=/data/GTEx/v8_jon/sample_batches/v8_batch00

#I usually run this script on a rugen node, but cd into tyrone so all output is kept there.


#Run a loop through the samples..
#Sorting the bam file
samtools sort -m 2G -o $file.readSorted.bam -O bam -n -@ 25 /data/GTEx/v8_jon/bams/$file.Aligned.sortedByCoord.out.patched.md.bam
# -m amount of memory per thread     -o output file     -O output file type     -n  sort by read name     -@ thread number 


#Converting bam file to fastq file
oedtools bamtofastq -i $file.readSorted.bam -fq $file.end1.fq -fq2 $file.end2.fq 2>temp_err_out.txt
#  -i input file    -fq1  output file for the first pair    -fq2  output file for the second pair
#This data is paired-end sequencing data, means there are 2 components for each read. Paired-end sequencing data is split into two files, one for each pair of the read
# I include 2>temp_err_out.txt because anytime there is an unpaired read, which can happen for various reasons, the tool will throw a warning, and any individual sample can have on the order of ~100,000s unpaired reads. 


#Mapping the sample and running the GATK pipeline produces a bunch of temporary and intermediate files that I just remove at the end
#Mapping the reads to the reference genome
STAR=/home/werner/tools/STAR-2.5.2b/source/STAR #This is the version of STAR Sara used previously
genome=/data/genomes/GRCh38_Gencode25  #Genome used for the already processed samples

read1=$file.end1.fq
read2=$file.end2.fq

$STAR --genomeDir $genome --readFilesIn $read1 $read2 --runThreadN 25 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM Unsorted --quantMode GeneCounts

#Rename the output file to indlude the sample name
mv Aligned.out.bam $file.Aligned.Unsorted.out.bam

samtools sort -m 2G -o $file.Aligned.SortedByCoord.bam -O bam  -@ 25 $file.Aligned.Unsorted.out.bam

# job 4. 1cpu, 

#Using the GATK pipeline to call variants in the sample
GATK=/home/werner/tools/gatk4.1.3.0/gatk/gatk
igvtools=/home/werner/tools/IGVTools/igvtools.jar
samtools=/home/werner/tools/samtools-1.9/samtools
genome_fa=/home/werner/xchrom_snp_skew/data/annot_genomes/GRCh38_Gencode25/GRCh38.p7.genome.fa #Genome fasta file location

chr='chrX'    #Specifies which chromosome to process, I'm only working with chrX 

chr_fa=/home/werner/xchrom_snp_skew/data/annot_genomes/GRCh38_Gencode25/$chr.fa   #Chromosome fasta file location

$samtools view -b $file.Aligned.sortedByCoord.out.bam $chr > $chr.bam    #Extract chromosome bam
$samtools view -b -q 10$chr.bam > $chr.filt.bam

echo "Adding read groups to bam file"
$GATK AddOrReplaceReadGroups \
    -I=$chr.filt.bam \
    -O=$chr.rg.bam \
    -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample

echo "Marking duplicates"
$GATK MarkDuplicates \
    -I=$chr.rg.bam  \
    -O=$chr.dedupped.bam \
    -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics

echo "Splitting and trimming"
$GATK SplitNCigarReads         \
   -R $genome_fa \
   -I $chr.dedupped.bam \
   -O $chr.split.filtered.bam


echo "Haplotype calling (?)"
$GATK HaplotypeCaller \
   -R $genome_fa \
   -L $chr \
   -I $chr.split.filtered.bam \
   --dont-use-soft-clipped-bases -stand-call-conf 20.0 \
   -O $chr.filtered.vcf

echo "Counting SNPs"
java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   $chr.split.filtered.bam   \
   $chr.split.filtered.wig   \
   $chr_fa

rm $chr.bam $chr.filt.bam $chr.rg.bam $chr.dedupped.bam $chr.dedupped.bai $chr.split.filtered.bam  $chr.split.filtered.bai output.metrics igv.log $file.Aligned.sortedByCoord.out.bam.bai

