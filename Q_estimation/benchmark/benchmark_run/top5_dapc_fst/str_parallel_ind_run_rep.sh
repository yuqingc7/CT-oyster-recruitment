#!/bin/bash -l

## start date
start=`date +%s`

mkdir -p logs output

## set environment
export PYTHONPATH=/programs/structure_threader-1.3.10/lib/python3.9/site-packages
export PATH=/programs/structure_threader-1.3.10/bin:$PATH

# Run STRUCTURE in parallel 4 more times
export STRUCTURE_MAINPARAMS="mainparams"
export STRUCTURE_EXTRAPARAMS="extraparams"

for rep in {1..4}; do
    echo "Starting replicate $rep..."
    mkdir -p logs/rep_${rep} output/rep_${rep}

    seed=$(( RANDOM + rep * 1000 ))
    echo "Using seed $seed for replicate $rep"

    parallel --jobs 24 /programs/structure_2_3_4/bin/structure \
        -m $STRUCTURE_MAINPARAMS -e $STRUCTURE_EXTRAPARAMS -D $seed -i {} \
        -o output/rep_${rep}/{/.} \> logs/rep_${rep}/{/.}.log 2>\&1 ::: structure_inputs/structure_run_*.strct_in

    echo "Finished replicate $rep."
done

# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$(((runtime % 3600) % 60))
echo "Total runtime: $hours:$minutes:$seconds (hh:mm:ss)"
