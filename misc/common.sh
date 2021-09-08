
nodeinfo() {
	echo "*********NODE*************************"
	hostname -f
	cat /etc/redhat-release
	NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
	echo "Processors: $NPROC "
	KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
	MBMEM=`expr $KMEM / 1000`
	echo "Memory MB: $MBMEM"
	NPROC=`cat /proc/cpuinfo  | grep processor | wc -l`
	echo "Processors: $NPROC "
	KMEM=`cat /proc/meminfo  | grep MemTotal | awk '{print $2}'`
	MBMEM=`expr $KMEM / 1000`
	echo "Memory MB: $MBMEM"
}

gettaskid(){
	if [ -z ${SGE_TASK_ID} ]; then 
		SGE_TASK_ID=1 
	fi
	echo $SGE_TASK_ID	
}

