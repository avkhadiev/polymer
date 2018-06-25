#!/bin/bash
nbonds=20
long=2
out=/Users/Arthur/stratt/polymer/test/gd_trial_job/out/
sim="${nbonds}_${long}"
declare -a names=(
                "${sim}_R" 
                "${sim}_rr"
                "${sim}_md" 
                "${sim}_sl_c"
                "${sim}_sl_l"
                "${sim}_p1_c"
                "${sim}_p1_l"
                "${sim}_p2_c"
                "${sim}_p2_l"
                )
out=/Users/Arthur/stratt/polymer/test/gd_trial_job/out/
for index in {0..8}
do
    nline=$(expr $index + 1)
    awk "FNR==${nline} { print; nextfile}" ${out}*out.log > ${names[index]}.txt
done
echo "aggregation complete"
