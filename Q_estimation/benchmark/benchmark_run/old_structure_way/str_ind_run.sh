#!/bin/bash -l

## start date
start=`date +%s`

## set environment
export PYTHONPATH=/programs/structure_threader-1.3.10/lib/python3.9/site-packages
export PATH=/programs/structure_threader-1.3.10/bin:$PATH

# Run STRUCTURE
/programs/structure_2_3_4/bin/structure \
-m mainparams_nocentroids -e extraparams -i test_selected5_905_nocentroids.recode.strct_in \
-o test_selected5_905_nocentroids &> test_selected5_905_nocentroids.log &

/programs/structure_2_3_4/bin/structure \
-m mainparams_onlycentroids -e extraparams -i test_selected5_905_onlycentroids.recode.strct_in \
-o test_selected5_905_onlycentroids &> test_selected5_905_onlycentroids.log &

/programs/structure_2_3_4/bin/structure \
-m mainparams_nocentroids_10k -e extraparams -i 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids.recode.strct_in \
-o 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids &> 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids.log &

/programs/structure_2_3_4/bin/structure \
-m mainparams_onlycentroids_10k -e extraparams -i 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids.recode.strct_in \
-o 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids &> 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids.log &

/programs/structure_2_3_4/bin/structure \
-m mainparams_nocentroids -e extraparams -i 905SNP.test_fp_noinvers.fixed.sorted_nocentroids.recode.strct_in \
-o 905SNP.test_fp_noinvers.fixed.sorted_nocentroids &> 905SNP.test_fp_noinvers.fixed.sorted_nocentroids.log &

/programs/structure_2_3_4/bin/structure \
-m mainparams_onlycentroids -e extraparams -i 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids.recode.strct_in \
-o 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids &> 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids.log &



# end date
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
