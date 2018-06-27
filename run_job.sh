#!/bin/bash
root=$HOME/polymer
task_id=$1
njobs=1000
chains=(2 20)
max_task_id=$(( ${#chains[@]} * $njobs ))
if [ "$task_id" -gt "$max_task_id" ]
then
	echo "$task_id is greaten than maximum ($max_task_id)"
else
	for idx in $(seq 1 ${#chains[@]})
	do
		istart=$(( ($idx - 1) * $njobs + 1))
		ifinish=$(( $istart + $njobs - 1 ))
		if [ "$task_id" -ge $istart -a "$task_id" -le $ifinish ]
		then
			nlinks=${chains[$(( $idx - 1 ))]}	
			cfg="${root}/configs/${nlinks}.cfg"
			sim=${nlinks}_$(( ((task_id - 1) % njobs) + 1 ))
			err="${root}/errs/${sim}.log"
			echo "Running MD"
			$root/bin/run_md_simulation $sim $cfg 1> /dev/null 2> $err 
			echo "MD ran, running geodesics" 
			log="${root}/logs/${sim}.log"
    			$root/bin/run_geodesic_simulation ${sim} ${cfg} 1> $log 2>$err
		fi
	done
fi
