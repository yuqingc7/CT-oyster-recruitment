# genotyped_merged_n1917_exmt.vcf
# Subsampled it to 1601

# excluding BAK23a HAM23a HOR23a LOP23a MCK23a MIH23a NRH23a 
# (1917-316=1601) 

vcftools --vcf genotyped_merged_n1917_exmt.vcf --keep sample_id_n1601_subset.txt --recode --recode-INFO-all --out genotyped_n1601_exmt
# After filtering, kept 1601 out of 1917 Individuals
# Outputting VCF file...
# After filtering, kept 65862 out of a possible 65862 Sites
# Run Time = 33.00 seconds

source /programs/miniconda3/bin/activate dDocent-2.8.13

./filter_missing_ind.sh genotyped_n1601_exmt.recode.vcf genotyped_n1601_exmt_indmiss
# Set cutoff 0.1, keep all individuals

# filter for 95% missingness
vcftools --vcf genotyped_n1601_exmt_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genotyped_n1601_exmt_maf05_maxmiss095
# After filtering, kept 1601 out of 1601 Individuals
# Outputting VCF file...
# After filtering, kept 54280 out of a possible 65862 Sites
# Run Time = 41.00 seconds

# Look for alleles that have very different missingness across batches
#Make a popmap file for two batches: cluster_n1601_batches.txt
awk '{count[$3]++} END {for (b in count) print b, count[b]}' cluster_n1601_batches.txt
#C 1244
#D 357
#Then separate list for batch C and D
vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf \
  --keep n1601_batch_C.txt \
  --missing-site \
  --out genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing

vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf \
  --keep n1601_batch_D.txt \
  --missing-site \
  --out genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing

# genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing.lmiss
# genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing.lmiss

#Extract SNP and missingness columns from both batches: 
awk 'NR>1 {print $1"_"$2 "\t" $6}' genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing.lmiss > genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing_miss.txt
awk 'NR>1 {print $1"_"$2 "\t" $6}' genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing.lmiss > genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing_miss.txt

join -1 1 -2 1 <(sort genotyped_n1601_exmt_maf05_maxmiss095_batch_C_missing_miss.txt) <(sort genotyped_n1601_exmt_maf05_maxmiss095_batch_D_missing_miss.txt) > joined_missingness.txt

# Calculate absolute difference and flag SNPs with big differences
awk '{diff = $2 - $3; if(diff < 0) diff = -diff; print diff}' joined_missingness.txt > joined_missingness_diff_values.txt

# Set 0.02 as the cutoff
awk -v c="$cutoff" '$2 > c {print $1}' joined_missingness_diff_with_ID.txt > SNPs_exceeding_cutoff.txt

wc -l SNPs_exceeding_cutoff.txt
# 4015 SNPs_exceeding_cutoff.txt

# Exclude SNPs from file
vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095.recode.vcf --exclude SNPs_exceeding_cutoff.txt \
--recode --recode-INFO-all --out genotyped_n1601_exmt_maf05_maxmiss095_missdiff"$cutoff"
# After filtering, kept 1601 out of 1601 Individuals
# Outputting VCF file...
# After filtering, kept 50265 out of a possible 54280 Sites
# Run Time = 36.00 seconds

# Need to make a popmap file
cut -f 2 popmap | sort | uniq | wc -l
# 42

./pop_missing_filter.sh genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02.recode.vcf popmap 0.95 42 genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095
# After filtering, kept 1601 out of 1601 Individuals
# Outputting VCF file...
# After filtering, kept 50265 out of a possible 50265 Sites
# Run Time = 34.00 seconds

# filter HWE by pop
./filter_hwe_by_pop.pl -v genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095.recode.vcf -p popmap -h 0.01 -c 0.5 -o genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe
# Outputting results of HWE test for filtered loci to 'filtered.hwe'
# Kept 48794 of a possible 50265 loci (filtered 1471 loci)

# Keep Biallelic only
/programs/plink-a2.3LM/plink2 --vcf genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe.recode.vcf --min-alleles 2 --max-alleles 2 --allow-extra-chr --make-bed --out genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic
/programs/plink-a2.3LM/plink2 --vcf genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe.recode.vcf --min-alleles 2 --max-alleles 2 --allow-extra-chr --export vcf --out genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic
# 48778 out of 48794 variants loaded from
# genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic-temporary.pvar.
# Note: No phenotype data present.
# 48778 variants remaining after main filters.

# LD clumping using r2 = 0.2
# Discarding 0 variant with MAC < 10 or MAF < 0.02.

# Run pruned.sh to generate pruned file: 
# After filtering, kept 1601 out of 1601 Individuals
# Outputting VCF file...
# After filtering, kept 38678 out of a possible 48794 Sites

# Remove inversions: 
vcftools --vcf genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic.pruned.recode.vcf --exclude-bed invers.bed --recode --recode-INFO-all --out genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic.pruned.nochr156invers
# After filtering, kept 38514 out of a possible 38678 Sites
