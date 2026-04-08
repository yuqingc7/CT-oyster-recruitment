!#/usr/bin/env bash

PREFIX="genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe"

vcftools --vcf $PREFIX".recode.vcf" --positions $PREFIX"_biallelic_LDclump_SNP.txt" --recode --recode-INFO-all --out $PREFIX"_biallelic.pruned"

/programs/plink2_linux_avx2_20230721/plink2  --vcf $PREFIX"_biallelic.pruned.recode.vcf" --make-bed --out $PREFIX"_biallelic.pruned" --allow-extra-chr

