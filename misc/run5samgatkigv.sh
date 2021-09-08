#!/bin/bash 
#   Usage: runsamgatvigv <setup> <infile> <outfile> <outvcf> <outwig> <outcut>
# 
#   args = $(setup) 
#          $(basefile)
#          $(chr)
#          $(genomedir)
#          $(outdir)
#					INFILE
#					/$(filebase).Aligned.sortedByCoord.out.bam  
#					OUTFILES
#          			/$(filebase).chrX.split.filtered.wig 
#                   /$(filebase).chrX.filtered.vcf 
#          			/$(filebase).chrX.intersect.wig 
#                   /$(filebase).chrX.cut
#   request_cpus = 1
#   request_memory = 10240
#
# job name
#$ -N run5samgatk
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-107 -tc 4
#
# processes per job
#$ -pe threads 1
#
#$ -wd /grid/gillis/data/hover/work/werner1/
#
# Per-thread memory request. 
#$ -l m_mem_free=10G
# 

#chr="chrX"
#chr_fa="$genomedir/$chr.fa"
#genome_fa="$genomedir/GRCh38.p7.genome.fa"

snpdb="/grid/gillis/home/hover/data/snpdb/chrX.snpdb.txt"

COMMON=~/git/cshl-work/scripts/werner1/common.sh 
echo "*********START*************************"
date
. $COMMON
nodeinfo


echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 5 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: run5samgatvigv <setup> <basefile> <workdir>"
    exit 1
fi

echo "Running job..."
echo "Running setup from $1"
. $1
IGV=`igvtools version`
echo "IGV version: $IGV"
echo "PATH=$PATH"

basefile=$2
# $SGE_TASK_ID
TASKID=$(gettaskid)
echo "TASKID is $TASKID"
filebase=`head -$TASKID $2 | tail -1 `
echo "Filebase is $filebase"

chr=$3
genomedir=$4
outdir=$5
 
# FILES
infile="$outdir/$filebase.Aligned.sortedByCoord.out.bam"
splitwig="$outdir/$filebase.$chr.split.filtered.wig" 
filtvcf="$outdir/$filebase.$chr.filtered.vcf" 
intwig="$outdir/$filebase.$chr.intersect.wig"
cut="$outdir/$filebase.$chr.cut"

mkdir -p $outdir/$filebase
cd $outdir/$filebase

echo "samtools view -b $infile $chr > $chr.bam"
samtools view -b $infile $chr > $chr.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "samtools view -b -q 10 $chr.bam > $chr.filt.bam"
samtools view -b -q 10 $chr.bam > $chr.filt.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo gatk AddOrReplaceReadGroups -I=$chr.filt.bam -O=$chr.rg.bam -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample
gatk AddOrReplaceReadGroups -I=$chr.filt.bam -O=$chr.rg.bam -SO=coordinate -RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Marking duplicates..."
echo gatk MarkDuplicates -I=$chr.rg.bam -O=$chr.dedupped.bam -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics
gatk MarkDuplicates -I=$chr.rg.bam -O=$chr.dedupped.bam -CREATE_INDEX=true -VALIDATION_STRINGENCY=SILENT -M=output.metrics
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Splitting and trimming..."
echo gatk SplitNCigarReads -R $genome_fa -I $chr.dedupped.bam -O $chr.split.filtered.bam
gatk SplitNCigarReads -R $genome_fa -I $chr.dedupped.bam -O $chr.split.filtered.bam
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Haplotype calling..."
echo gatk HaplotypeCaller -R $genome_fa -L $chr -I $chr.split.filtered.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0  -O $chr.filtered.vcf
gatk HaplotypeCaller -R $genome_fa -L $chr -I $chr.split.filtered.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0  -O $chr.filtered.vcf
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

TMPFILE=`mktemp -p ./  tmp.XXXXXX.wig`
# #chr_fa="$genomedir/$chr.fa"
echo "TMPFILE is $TMPFILE"
echo "Counting SNPs..."
echo igvtools count -z 0 -w 1 --bases --strands read $chr.split.filtered.bam $TMPFILE $chr_fa
igvtools count -z 0 -w 1 --bases --strands read $chr.split.filtered.bam $TMPFILE $chr_fa
RET=$?
echo "Command Return code was $RET"
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Additional filtering/cutting..."
grep $chr $chr.filtered.vcf | cut -f2-5 > A
cut -f1 A > B
grep -Fw -f B $snpdb > D
grep -Fw -f D $TMPFILE > $chr.intersect.wig
grep -Fw -f D A > $chr.cut 


echo "**************FINALIZE************"
ls -alh

# handle TMPFILE (which is chrX.split.filtered.wig) 
echo mv -v $TMPFILE $splitwig
time mv -v $TMPFILE $splitwig

# handle files...
echo mv -v $chr.filtered.vcf $filtvcf
mv -v $chr.filtered.vcf $filtvcf

echo mv -v $chr.intersect.wig $intwig
mv -v $chr.intersect.wig $intwig

echo mv -v $chr.cut $cut
mv -v $chr.cut $cut

cd ..
rm -rf $filebase/*
rmdir $filebase 

echo "*********END***************************"
date
sleep 30
exit $RET

