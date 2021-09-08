#!/bin/bash 
#   Usage: runsamtools <setup> <indir> <outdir>
#
# args = $(setup)
#        $(indir)      # in:   /$(filebase).Aligned.sortedByCoord.out.patched.md.bam 
#        $(outdir)     # out:  /$(filebase).readSorted.bam 
# request_cpus = 12
# request_memory = 16348
##
# job name
#$ -N run1samtools
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-357 -tc 10
#
# processes per job
#$ -pe threads 8 
#
#$ -wd /grid/gillis/data/hover/work/werner1/
#
# Per-thread memory request. 
#$ -l m_mem_free=3G
# 


echo "*********START*************************"
date

echo "*********NODE*************************"
hostname -f
cat /etc/redhat-release

NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
echo "Processors: $NPROC "
KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
MBMEM=`expr $KMEM / 1000`
echo "Memory MB: $MBMEM"

echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 4 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runsamtools <setup> <basemap> <infile> <outfile>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

basefile=$2
# $SGE_TASK_ID
filebase=`head -$SGE_TASK_ID $2 | tail -1 `
echo "Filebase is $filebase"
infile="$3/$filebase.Aligned.sortedByCoord.out.patched.md.bam"
outfile="$4/$filebase.readSorted.bam"


echo "Running job..."
echo samtools sort -m 2G -o $outfile -O bam -n  -@ 8 $infile
time samtools sort -m 2G -o $outfile -O bam -n -@ 8 $infile
RET=$?
echo "Job command Return code was $RET"

echo "*********DONE*************************"
date
echo "*********END***************************"
exit $RET

