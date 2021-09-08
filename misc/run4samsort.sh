#!/bin/bash 
#   Usage: runsamtools <setup> <infile> <outfile>
#
#   args = $(setup) 
#          $(basefile)
#          $(outdir)

#          			$(outdir)/$(filebase).Aligned.out.bam  
#          			$(outdir)/$(filebase).Aligned.sortedByCoord.out.bam
#   request_cpus = 10
#   request_memory = 20480
#
# job name
#$ -N run4sam
#
# job indexes for array all jobs, but only run 10 at a time for disk quota 
#$ -t 1-357 -tc 4
#
# processes per job
#$ -pe threads 20
#
#$ -wd /grid/gillis/data/hover/work/werner1/
#
# Per-thread memory request. 
#$ -l m_mem_free=2G
# 
COMMON=~/git/cshl-work/scripts/werner1/common.sh 
NUMTHR=20

echo "*********START*************************"
date
. $COMMON
nodeinfo

echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 3 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: run4samsort <setup> <basefile> <workdir>"
    exit 1
fi

echo "Running setup from $1"
. $1
echo "PATH=$PATH"

basefile=$2
# $SGE_TASK_ID
TASKID=$(gettaskid)
echo "TASKID is $TASKID"
filebase=`head -$TASKID $2 | tail -1 `
echo "Filebase is $filebase"

outdir=$3
infile="$outdir/$filebase.Aligned.out.bam"
outfile="$outdir/$filebase.Aligned.sortedByCoord.out.bam"

echo "Running job..."
#TMPFILE=`mktemp -p ./`
#echo "TMPFILE is $TMPFILE"
echo samtools sort -m 2G -o $outfile -O bam  -@ $NUMTHR $infile
time samtools sort -m 2G -o $outfile -O bam  -@ $NUMTHR $infile
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo samtools index $outfile 
time samtools index $outfile 
RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi

echo "Job command Return code was $RET"

echo "*********DONE*************************"
date

echo "*********END***************************"
exit $RET

