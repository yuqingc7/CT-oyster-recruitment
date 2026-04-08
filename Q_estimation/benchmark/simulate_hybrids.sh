python recom-sim_yc.py --num-offs 20 --p1name Wild --p2name AQ --exclude --out in_silico_hybrid3_fp_noinvers geno_HR_AQ_n284_fp_noinvers.genepop2 3
# This gives us - F1HYB, B2Wild, B3Wild

# Now manually run Wild x B3Wild to create B4Wild
python recom-sim_yc.py --num-offs 20 --p1name Wild --p2name B3Wild --exclude --out in_silico_B4Wild_fp_noinvers in_silico_HRxB3Wild_fp_noinvers.gen 1

# Manually run Wild x B4Wild to create B5Wild
python recom-sim_yc.py --num-offs 20 --p1name Wild --p2name B4Wild --exclude --out in_silico_B5Wild_fp_noinvers in_silico_HRxB4Wild_fp_noinvers.gen 1

# Manually run Wild x Wild to create WildUnadmixed
python recom-sim_yc.py --num-offs 20 --p1name Wild --p2name Wild --exclude --out in_silico_WildUnadmixed_fp_noinvers geno_HR_HR_n50x2_fp_noinvers.genepop2 1

# We need a test gen file containing: 
# F1HYB = 50% 
# 	n=20, in_silico_hybrid3_fp_noinvers.gen
# B2Wild = 25%
# 	n=20, in_silico_hybrid3_fp_noinvers.gen
# B3Wild = 12.5%
# 	n=20, in_silico_HRxB3Wild_fp_noinvers.gen
# B4POP1 = manually run Wild x B3POP1 = 6.25%
# 	n=20, in_silico_B4Wild_fp_noinvers.gen
# B5POP1 = manually run Wild x B4POP1 = 3.125%
# 	n=20, in_silico_B5Wild_fp_noinvers.gen
# Unadmixed = manually run Wild x Wild = 0%
# n=20, in_silico_WildUnadmixed_fp_noinvers.gen
# AQ = 100%
# n=234, geno_HR_AQ_n284_fp_noinvers.genepop2
# HR = 0%
# n=50, geno_HR_AQ_n284_fp_noinvers.genepop2

# Centroid wild (Wild x Wild) = 0% -> Ref
# n=100, in_silico_wild_AQ_fp_noinvers.gen
# Centroid AQ (AQ x AQ) =100% -> Ref 
# n=100, in_silico_wild_AQ_fp_noinvers.gen

# Total 604 individuals 

sed -i 's/ ,  /,  /g' test_fp_noinvers.gen
tail -n 202 in_silico_wild_AQ_fp_noinvers.gen >> test_fp_noinvers.gen
tail -n 21 in_silico_WildUnadmixed_fp_noinvers.gen >> test_fp_noinvers.gen
tail -n 20 in_silico_B5Wild_fp_noinvers.gen >> test_fp_noinvers.gen
tail -n 20 in_silico_B4Wild_fp_noinvers.gen >> test_fp_noinvers.gen
tail -n 20 in_silico_HRxB3Wild_fp_noinvers.gen >> test_fp_noinvers.gen
grep "B2Wild_" in_silico_hybrid3_fp_noinvers.gen >> test_fp_noinvers.gen
grep "F1HYB_" in_silico_hybrid3_fp_noinvers.gen >> test_fp_noinvers.gen

# https://thierrygosselin.github.io/radiator/reference/genomic_converter.html 
# Convert to VCF formats
# Run convert.R
# -> test_fp_noinvers.vcf

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
}' test_fp_noinvers.vcf > test_fp_noinvers.fixed.vcf
