# Subset loci
bcftools view -i 'ID=@../SNP_selected_5.list' test_fp_noinvers.fixed.sorted.vcf -o test_selected5_905.vcf

# Convert to structure format
/programs/plink-1.9-x86_64-beta5/plink --vcf test_selected5_905.vcf \
--allow-extra-chr --make-bed \
--out test_selected5_905

/programs/plink-1.9-x86_64-beta5/plink --bfile test_selected5_905 \
--allow-extra-chr \
--recode structure --out /local/workdir/yc2644/CV_CT_array/structure_test/test_recomsim/test_selected5_905

# Test on Simulated data
# Run all inds together with centroids refs, 905 markers
/programs/structure_2_3_4/bin/structure \
> -m mainparams_onlycentroids -e extraparams -i test_selected5_905_onlycentroids.recode.strct_in \
> -o test_selected5_905_onlycentroids &> test_selected5_905_onlycentroids.log &

# Run all inds together with empirical refs, 905 markers
/programs/structure_2_3_4/bin/structure \
> -m mainparams_nocentroids -e extraparams -i test_selected5_905_nocentroids.recode.strct_in \
> -o test_selected5_905_nocentroids &> test_selected5_905_nocentroids.log &

# Run all inds together with empirical refs, 10k markers
# randomly subsample 10k SNPs from vcf files
bcftools view --header-only  test_fp_noinvers.fixed.sorted.vcf > 10kSNP.test_fp_noinvers.fixed.sorted.vcf
bcftools view --no-header  test_fp_noinvers.fixed.sorted.vcf | awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t' -k1,1g | cut -f2-  | head -n 10000 >> 10kSNP.test_fp_noinvers.fixed.sorted.vcf

/programs/plink-1.9-x86_64-beta5/plink --vcf 10kSNP.test_fp_noinvers.fixed.sorted.vcf \
--allow-extra-chr --make-bed \
--out 10kSNP.test_fp_noinvers.fixed.sorted

/programs/plink-1.9-x86_64-beta5/plink --bfile 10kSNP.test_fp_noinvers.fixed.sorted \
--allow-extra-chr \
--recode structure --out /local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/old_structure_way/10kSNP.test_fp_noinvers.fixed.sorted

/programs/structure_2_3_4/bin/structure \
> -m mainparams_nocentroids_10k -e extraparams -i 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids.recode.strct_in \
> -o 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids &> 10kSNP.test_fp_noinvers.fixed.sorted_nocentroids.log &

# Run all inds together with centroid refs, 10k markers
/programs/structure_2_3_4/bin/structure \
> -m mainparams_onlycentroids_10k -e extraparams -i 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids.recode.strct_in \
> -o 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids &> 10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids.log &

# Run all inds together with empirical refs, 905 random markers
# randomly subsample 905 SNPs from vcf files
bcftools view --header-only  test_fp_noinvers.fixed.sorted.vcf > 905SNP.test_fp_noinvers.fixed.sorted.vcf
bcftools view --no-header  test_fp_noinvers.fixed.sorted.vcf | awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t' -k1,1g | cut -f2-  | head -n 905 >> 905SNP.test_fp_noinvers.fixed.sorted.vcf

/programs/plink-1.9-x86_64-beta5/plink --vcf 905SNP.test_fp_noinvers.fixed.sorted.vcf \
--allow-extra-chr --make-bed \
--out 905SNP.test_fp_noinvers.fixed.sorted

/programs/plink-1.9-x86_64-beta5/plink --bfile 905SNP.test_fp_noinvers.fixed.sorted \
--allow-extra-chr \
--recode structure --out /local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/old_structure_way/905SNP.test_fp_noinvers.fixed.sorted

/programs/structure_2_3_4/bin/structure \
> -m mainparams_nocentroids -e extraparams -i 905SNP.test_fp_noinvers.fixed.sorted_nocentroids.recode.strct_in \
> -o 905SNP.test_fp_noinvers.fixed.sorted_nocentroids &> 905SNP.test_fp_noinvers.fixed.sorted_nocentroids.log &

# Run all inds together with centroid refs, 905 random markers
/programs/structure_2_3_4/bin/structure \
> -m mainparams_onlycentroids -e extraparams -i 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids.recode.strct_in \
> -o 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids &> 905SNP.test_fp_noinvers.fixed.sorted_onlycentroids.log &

# Run one by one with centroid refs, 905 random markers
