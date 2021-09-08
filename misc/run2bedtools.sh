#!/bin/bash 
#   Usage: runbedtools <setup> <infile> <outfile1> <outfile2>
#
#   args =  $(setup)
#           $(basefile)	
#           $(outdir)				/$(filebase).readSorted.bam 
#           						$(outdir)/$(filebase).end1.fq  
#           						$(outdir)/$(filebase).end2.fq
#   request_cpus = 1
#   request_memory = 2048

# job name
#$ -N run2bedtools
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-357 -tc 10
#
# processes per job
#$ -pe threads 1 
#
#$ -wd /grid/gillis/data/hover/work/werner1/
#
# Per-thread memory request. 
#$ -l m_mem_free=2G
# 
COMMON=~/git/cshl-work/scripts/werner1/common.sh 

echo "*********START*************************"
date
. $COMMON
echo "*********JOB*************************"

echo "Args are $@"
if [ $# -ne 3 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runbedtools <setup> <infile> <outfile1> <outfile2>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

echo "Running job..."

basefile=$2
# $SGE_TASK_ID
filebase=`head -$SGE_TASK_ID $2 | tail -1 `
echo "Filebase is $filebase"
infile="$3/$filebase.readSorted.bam"
outfile1="$3/$filebase.end1.fq "
outfile2="$3/$filebase.end2.fq"

echo bedtools bamtofastq -i $infile -fq $outfile1 -fq2 $outfile2 
time bedtools bamtofastq -i $infile -fq $outfile1 -fq2 $outfile2
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "*********DONE*************************"
date
echo "*********END***************************"

exit $RET
