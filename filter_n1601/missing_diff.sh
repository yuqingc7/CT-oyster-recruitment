#!/bin/bash

# Make a popmap file for two batches: cluster_n1601_batches.txt
# awk '{count[$3]++} END {for (b in count) print b, count[b]}' cluster_n1601_batches.txt
# C 1244
# D 357
# Then separate list for batch C and D

vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf \
  --keep n1601_batch_C.txt \
  --missing-site \
  --out genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing

vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf \
  --keep n1601_batch_D.txt \
  --missing-site \
  --out genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing

# Extract SNP and missingness columns from both batches: 
awk 'NR>1 {print $1"_"$2 "\t" $6}' genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing.lmiss > genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing_miss.txt
awk 'NR>1 {print $1"_"$2 "\t" $6}' genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing.lmiss > genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing_miss.txt

join -1 1 -2 1 <(sort genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing_miss.txt) <(sort genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing_miss.txt) > joined_missingness.txt

# Calculate absolute difference
awk '{diff = $2 - $3; if(diff < 0) diff = -diff; print $1, diff}' joined_missingness.txt > joined_missingness_diff_with_ID.txt
awk '{diff = $2 - $3; if(diff < 0) diff = -diff; print diff}' joined_missingness.txt > joined_missingness_diff_values.txt

# Plot the histogram of absolute differences
# Using gnuplot to visualize the differences
gnuplot << \EOF
set terminal dumb size 120, 30
set title "Histogram of Missingness Differences Between Batches"
set xlabel "Absolute Missingness Difference"
set ylabel "Number of SNPs"
set lmargin 10
set autoscale

binwidth=0.001
bin(x,width)=width*floor(x/width) + width/2.0

plot 'joined_missingness_diff_values.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes

pause -1
\EOF

# Pick a cutoff for the absolute difference from the histogram
echo "Please enter a cutoff for the absolute difference in missingness between batches (e.g., 0.05):"
read cutoff

#  and flag SNPs with big differences
awk -v c="$cutoff" '$2 > c {print $1}' joined_missingness_diff_with_ID.txt > SNPs_exceeding_cutoff.txt

# Create a list of SNPs with significant differences
vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf --exclude SNPs_exceeding_cutoff.txt \
--recode --recode-INFO-all --out genotyped_n1601_exmt_maf05_maxmiss095_missdiff"$cutoff"
