#!/bin/bash -l

## start date
start=`date +%s`

#PREFIX=10kSNP.genotyped_n1917_fp_noinvers
SAMPLE=n1601

WORKDIR=/workdir/yc2644/CV_CT_array

## set environment
export PYTHONPATH=/programs/structure_threader-1.3.10/lib/python3.9/site-packages
export PATH=/programs/structure_threader-1.3.10/bin:$PATH

mkdir $WORKDIR"/structure/"$SAMPLE"_cor"

# run software
structure_threader run -i $WORKDIR"/structure/threader_input/"$SAMPLE".recode.strct_in" \
-o $WORKDIR"/structure/"$SAMPLE"_cor" -st /programs/structure_2_3_4/bin/structure -Klist 2 3 4 5 6 7 8 9 10 -R 10 -t 24 --log true


# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
