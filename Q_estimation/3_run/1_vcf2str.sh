# Subset to informative loci
bcftools view -i 'ID=@../SNP_selected_5.list' geno_n1601_in_silico_wild_AQ_fp_noinvers.fixed.vcf -o geno_n1601_in_silico_wild_AQ_selected5_905.vcf

# Convert to structure format
/programs/plink-1.9-x86_64-beta5/plink --vcf geno_n1601_in_silico_wild_AQ_selected5_905.vcf \
--allow-extra-chr --make-bed \
--out geno_n1601_in_silico_wild_AQ_selected5_905

/programs/plink-1.9-x86_64-beta5/plink --bfile geno_n1601_in_silico_wild_AQ_selected5_905 \
--allow-extra-chr \
--recode structure --out /local/workdir/yc2644/CV_CT_array/structure_K/geno_n1601_in_silico_wild_AQ_selected5_905

# Run one at a time with centroid refs
bash str_prepare_ind_input.sh

bash str_parallel_ind_run.sh &
