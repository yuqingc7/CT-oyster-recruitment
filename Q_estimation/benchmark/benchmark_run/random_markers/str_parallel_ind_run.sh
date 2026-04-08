#!/bin/bash -l

## start date
start=`date +%s`

## set environment
export PYTHONPATH=/programs/structure_threader-1.3.10/lib/python3.9/site-packages
export PATH=/programs/structure_threader-1.3.10/bin:$PATH

# Run STRUCTURE in parallel
export STRUCTURE_MAINPARAMS="mainparams"
export STRUCTURE_EXTRAPARAMS="extraparams"

mkdir -p output
mkdir -p logs

parallel --jobs 24 /programs/structure_2_3_4/bin/structure \
-m $STRUCTURE_MAINPARAMS -e $STRUCTURE_EXTRAPARAMS -i {} \
-o output/{/.} '&> logs/{/.}.log' ::: structure_inputs/structure_run_*.strct_in

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
