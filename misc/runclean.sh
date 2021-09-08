#!/bin/bash 
#   Usage: runclean <setup> <outdir> <filebase> 
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
CV=`condor_version | tr -d "\n"`
echo "Condor Version: $CV"
NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
echo "Processors: $NPROC "
KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
MBMEM=`expr $KMEM / 1000`
echo "Memory MB: $MBMEM"

echo "*********JOB*************************"
echo "Args are $@"
if [ $# -ne 3 ]; then
    echo "Incorrect number of arguments."
    echo "Usage: runclean <setup> <outdir> <filebase>" 
    exit 1
fi

echo "Running setup from $1"
. $1
echo $PATH
echo "Running job..."

# GTEX-1J8JJ-0126-SM-D3L93.chrX.filtered.vcf
# GTEX-1J8JJ-0126-SM-D3L93.chrX.split.filtered.wig
# GTEX-1J8JJ-0126-SM-D3L93.chrX.intersect.wig
# GTEX-1J8JJ-0126-SM-D3L93.chrX.cut

# To delete:
#<base>.end1.fq
#<base>.end2.fq
#<base>.readSorted.bam
#<base>.Aligned.out.bam

if [ -s $2/$3.chrX.filtered.vcf ] && [ -s $2/$3.chrX.split.filtered.wig ] && [ -s  $2/$3.chrX.intersect.wig ] && [ -s  $2/$3.chrX.cut ]; then
    echo "Processing successful, cleaning up..."
    ENDS=" end1.fq  end2.fq readSorted.bam Aligned.out.bam"
    for end in $ENDS; do
        echo rm -f $2/$3.$end
        rm -f $2/$3.$end
    done	
else
    echo "Required files empty or missing. Problem..."
fi	

RET=$?
if [ $RET -ne 0 ] ; then
    exit $RET
fi
echo "Job command Return code was $RET"

echo "*********DONE*************************"
date
echo "*********END***************************"

exit $RET
