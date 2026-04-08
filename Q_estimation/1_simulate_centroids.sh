# Simulate centroids for in silico crosses of HR and AQ, using the same 38,514 SNPs as in the real data.
python recom-sim_yc.py --num-offs 100 --p1name Wild --p2name Wild --exclude --out in_silico_wild_fp_noinvers geno_HR_HR_n50x2_fp_noinvers.genepop2 1
python recom-sim_yc.py --num-offs 100 --p1name AQ --p2name AQ --exclude --out in_silico_AQ_fp_noinvers geno_AQ_AQ_n234x2_fp_noinvers.genepop2 1

# rename files for later use
mv *txt *gen 

# Replace every individual name F1HYB_ -> WILD;  F1HYB_ -> AQ
sed -i 's/F1HYB_/WILD_/g' in_silico_wild_fp_noinvers.gen
sed -i 's/F1HYB_/AQ_/g' in_silico_AQ_fp_noinvers.gen

# Put Gen files together into one
# -> in_silico_wild_AQ_fp_noinvers.gen

# Now turn to the real data gen file, which contains 1601 individuals (38,514 SNPs)
# Unify “,  “
sed 's/, /,  /g' geno_n1601_fp_noinvers.genepop2 > geno_n1601_fp_noinvers.genepop3
# Also replace - with _
sed -i 's/-/_/g' geno_n1601_fp_noinvers.genepop3

# Combine in_silico_wild_AQ_fp_noinvers.gen and geno_n1601_fp_noinvers.genepop3
# -> geno_n1601_in_silico_wild_AQ_fp_noinvers.gen

# https://thierrygosselin.github.io/radiator/reference/genomic_converter.html 
# Convert to VCF formats
# Run convert.R
# -> geno_n1601_in_silico_wild_AQ_fp_noinvers.vcf

# Keep all ## header lines and #CHROM header intact,
# CHROM = the first number before _ in the ID (10),
# POS = the second number (1000149),
# ID = 10_1000149
awk 'BEGIN{FS=OFS="\t"} 
/^##/ {print; next} 
/^#CHROM/ {print; next} 
{
  split($3, a, "_"); 
  $1 = a[1]; 
  $2 = a[2]; 
  $3 = a[1] "_" a[2]; 
  print
}' geno_n1601_in_silico_wild_AQ_fp_noinvers.vcf > geno_n1601_in_silico_wild_AQ_fp_noinvers.fixed.vcf

# Add contigs to the header
    ##contig=<ID=10>
    ##contig=<ID=1>
    ##contig=<ID=2>
    ##contig=<ID=3>
    ##contig=<ID=4>
    ##contig=<ID=5>
    ##contig=<ID=6>
    ##contig=<ID=7>
    ##contig=<ID=8>
    ##contig=<ID=9>

