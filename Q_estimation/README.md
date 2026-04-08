# Estimate introgression levels

1. Simulate in silico centroid references using recom-sim
   - Optional: Evaluate your centroids
2. Identify informative loci and subset vcf to those loci
3. Run one individual at a time with centroid references
   - Prepare input (`1_vcf2str.sh`, `2_str_prepare_ind_input.sh`)
   - Obtain Pind from structure runs (`3_str_parallel_ind_run.sh`)
   - Obtain introgression levels by normalizing Pind using P_AQ and P_Wild (averages of empirical reference samples) (`4_sum_Q_test.R`)


## Optional: Benchmark your approah (`benchmark` folder)

1. Simulate a range of hybrid classes
2. Apply your approach (see `benchmark_approaches.sh` for examples)
3. Compare your estimated introgression levels against the expectations (see `sum_Q_test.R` and `sum_Q_test_old_structure.R`)
