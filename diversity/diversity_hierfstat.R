setwd("/local/workdir/yc2644/CV_CT_array/hierfstat")

library(tidyverse)
library(MASS)
#install_github("jgx65/hierfstat")
library("hierfstat")
library(gaston,quietly=TRUE)
library("adegenet")
library(vcfR)

# # import VCF files and convert to genind
# filename <- "10kSNP.geno_n1601_f.vcf"
# filepath <- paste0("/local/workdir/yc2644/CV_CT_array/vcf_n1601/", filename)
# all_vcf <- read.vcfR(filepath, verbose = FALSE)
# all_genind <- vcfR2genind(all_vcf)

filename <- "10kSNP.geno_n1601_fp_noinvers.vcf"
filepath <- paste0("/local/workdir/yc2644/CV_CT_array/vcf_n1601/", filename)
pruned_vcf <- read.vcfR(filepath, verbose = FALSE)
pruned_genind <- vcfR2genind(pruned_vcf)

# load population information
pop_info <- read_tsv("../vcf_n1601/pop_grp_info_n1601.tsv", col_names = T)
pop_info$Pop = factor(pop_info$Pop)
pop_info$Location = factor(pop_info$Location)
pop_info$Grp = factor(pop_info$Grp)
pop_info$Type = factor(pop_info$Type)


# basic stats -------------------------------------------------------------

all_genind@pop <- pop_info$Pop
#all_genind@pop <- pop_info$locat
all_genind

# `basic.stats` (it calculates $H_O$, $H_S$, $F_{IS}$, $F_{ST}$ etc...)
basicstat <- basic.stats(all_genind, diploid = TRUE, digits = 2)
names(basicstat)
basicstat$overall

# mean Ho per population
Ho <- colMeans(x=basicstat$Ho, na.rm = TRUE)
# wilcox.test(c(Ho[1],Ho[4:19]),Ho[2:3], alternative = "two.sided")
write.table(Ho, file = "10k.n1601.pop.Ho.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# mean He per population
He <- colMeans(x=basicstat$Hs, na.rm = TRUE)
write.table(He, file = "10k.n1601.pop.He.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# mean Fis per population
Fis <- colMeans(x=basicstat$Fis, na.rm = TRUE)
Fis_ci <- boot.ppfis(all_genind)$fis.ci
cbind(Fis, Fis_ci)
write.table(cbind(Fis, Fis_ci), file = "10k.n1601.pop.Fis.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# allelic richness
Arich <- allelic.richness(all_genind,min.n=NULL,diploid=TRUE)
ind_mean <- colMeans(x=Arich$Ar, na.rm = TRUE)
ind_mean
write.table(ind_mean, file = "10k.n1601.individual.allelic.richness.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)


# Functions beta.dosage give estimates of individual inbreeding coefficients
# and kinships between individuals
# convert to a dosage format via fstat2dos
all_genind.dos <- fstat2dos(all_genind) #converts fstat format to dosage
fis.all <- fis.dosage(all_genind.dos, pop=all_genind@pop)
beta.dos.all <- beta.dosage(all_genind.dos)
# By default (inb=TRUE) the inbreeding coefficient
# is returned on the main diagonal.
indF.all <- pop_info %>% mutate(indF = diag(beta.dos.all))
write.table(indF.all, file = "10k.n1601.indF.all.txt", sep = "\t", quote = FALSE,
            row.names = F, col.names = T)

indF.pop <- indF.all %>%
  group_by(Pop) %>%
  summarise(average=mean(indF))
write.table(indF.pop, file = "10k.n1601.indF.pop.txt", sep = "\t", quote = FALSE,
            row.names = F, col.names = T)

# Fst -------------------------------------------------------------

pruned_genind@pop <- pop_info$Pop
pruned_genind
# Weir and Cockerham's estimate
#wc(pruned_genind)

# Pairwise Fst
fst_pruned <- genet.dist(pruned_genind, method = "WC84")
# pruned_df <- genind2hierfstat(pruned_genind)
# fst_pruned <- pairwise.WCfst(pruned_df)

write.matrix(fst_pruned, file = "10k.n1601.pairwise.pruned.fst.txt", sep = "\t")

