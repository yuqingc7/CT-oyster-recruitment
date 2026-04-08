# /local/workdir/yc2644/CV_CT_array/lists/HR_n50.txt
# /local/workdir/yc2644/CV_CT_array/lists/AQ_n234.txt
# /local/workdir/yc2644/CV_CT_array/lists/CT_ER_new_n1317.txt

# geno_n1601_fp_noinvers.vcf (38,514 SNPs)

bcftools view -S /local/workdir/yc2644/CV_CT_array/lists/HR_n50.txt -Ov -o geno_HR_n50_fp_noinvers.vcf geno_n1601_fp_noinvers.vcf
bcftools view -S /local/workdir/yc2644/CV_CT_array/lists/AQ_n234.txt -Ov -o geno_AQ_n234_fp_noinvers.vcf geno_n1601_fp_noinvers.vcf
cat /local/workdir/yc2644/CV_CT_array/lists/HR_n50.txt /local/workdir/yc2644/CV_CT_array/lists/AQ_n234.txt > /local/workdir/yc2644/CV_CT_array/lists/HR_AQ_n284.txt
bcftools view -S /local/workdir/yc2644/CV_CT_array/lists/HR_AQ_n284.txt -Ov -o geno_HR_AQ_n284_fp_noinvers.vcf geno_n1601_fp_noinvers.vcf

# install vcf2popgen and dependencies
conda create -n vcf2popgen python=3.9
conda activate vcf2popgen
pip install vcf2popgen-0.2.1.tar.gz
pip install "numpy<1.27" matplotlib pandas scipy pandas-plink jsonschema
python -m pip check

# convert to genepop format (in python)
python3
import vcf2popgen
data = vcf2popgen.read('geno_HR_AQ_n284_fp_noinvers.vcf', '../lists/HR_AQ_n284_map.txt')
data.to_genepop('geno_HR_AQ_n284_fp_noinvers.genepop')

data = vcf2popgen.read('geno_n1601_fp_noinvers.vcf', '../lists/ALL_n1601_map.txt')
data.to_genepop('geno_n1601_fp_noinvers.genepop')
exit

# remove snp_ from each loci lines
sed -E 's/snp[0-9]+_//g' geno_HR_AQ_n284_fp_noinvers.genepop > geno_HR_AQ_n284_fp_noinvers.genepop2
sed -E 's/snp[0-9]+_//g' geno_n1601_fp_noinvers.genepop > geno_n1601_fp_noinvers.genepop2
