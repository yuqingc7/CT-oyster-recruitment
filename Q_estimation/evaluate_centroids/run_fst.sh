# Define populations
#pops=("geno_AQ_samples.txt" "geno_HR_samples.txt" "geno_in_silico_AQ_samples.txt" "geno_in_silico_wild_samples.txt")
pops=(geno_AQ1_samples.txt geno_AQ2_samples.txt geno_AQ3_samples.txt geno_AQ4_samples.txt geno_AQ5_samples.txt geno_AQ6_samples.txt geno_HR_samples.txt geno_in_silico_AQ_samples.txt geno_in_silico_wild_samples.txt)

# Loop over all pairs
n=${#pops[@]}
for ((i=0; i<n-1; i++)); do
  for ((j=i+1; j<n; j++)); do
    pop1=${pops[$i]}
    pop2=${pops[$j]}
    outname=$(basename $pop1 .txt)_vs_$(basename $pop2 .txt)
    vcftools --vcf geno_n1601_in_silico_wild_AQ_selected5_905.vcf \
             --weir-fst-pop $pop1 \
             --weir-fst-pop $pop2 \
             --out $outname
  done
done

grep "weighted" *_vs_*.log > mean_FST_all.txt
