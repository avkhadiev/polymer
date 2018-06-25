#!/bin/bash
nbonds=20
long=2
cfg=/Users/Arthur/stratt/polymer/configs/trial_job.cfg
out=/Users/Arthur/stratt/polymer/test/gd_trial_job/out/
for njob in {1..100}
do
    sim="${nbonds}_${long}_${njob}"
    log1="${out}/${sim}_gd_out.log"
    log2="${out}/${sim}_gd_err.log"
    echo "Running job ${njob}"
    bin/run_geodesic_simulation ${sim} ${cfg} 1>${log1} 2>${log2}
done

