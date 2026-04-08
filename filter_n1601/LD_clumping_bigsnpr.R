#setwd("/local/workdir/yc2644/CV_CT_array/check_inversion/filtering")

setwd("/local/workdir/yc2644/CV_CT_array/refiltering_n1601")

library(bigsnpr)

prefix="genotyped_n1601_exmt_maf05_maxmiss095_missdiff0.02_popmiss095_hwe_biallelic"
NCORES <- nb_cores()

snp_readBed2(paste(prefix,".bed", sep=""), 
                        ncores = NCORES)
obj.bigSNP <- snp_attach(paste(prefix,".rds", sep=""))

G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

# check if there is any missing values as NA
big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
# genotype imputation
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
big_counts(G, ind.col = 1:dim(G)[1]) # check if NAs are 0

# LD clumping using r2 = 0.2
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, 
                     ncores = NCORES, thr.r2 = 0.2, size = 10) # size is the window size of 10K

# extract SNPs after clumping
which_pruned = attr(newpc, which="subset")
#snp_list <- SNPs[which_pruned]
chr_list <- CHR[which_pruned]
pos_list <- POS[which_pruned]

print(paste0("SNPs after clumping is: ", length(pos_list), " out of ", dim(obj.bigSNP$map)[1]))

# generate LD clumped SNP list
snp_list <- tibble(chr_list, pos_list)
write.table(snp_list, file = paste0(prefix, "_LDclump_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# run pruned_vcf.sh to generate pruned VCF file & BED file
