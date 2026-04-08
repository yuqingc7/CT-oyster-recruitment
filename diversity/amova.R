setwd("/local/workdir/yc2644/CV_CT_array/hierfstat")

library(tidyverse)
library(MASS)
#install_github("jgx65/hierfstat")
library("hierfstat")
library(gaston,quietly=TRUE)
library("adegenet")
library(vcfR)

# AMOVA -------------------------------------------------------------
library(poppr)

# strata(pruned_genind) <- pop_info
# pruned_genind
# amova.result <- poppr.amova(pruned_genind, ~grp/locat, nperm=999)
# amova.test <- randtest(amova.result,nrepet = 999) # Test for significance

# subset to Adult vs Spat contrast
filename <- "10kSNP.geno_n1601_fp_noinvers_contrast.recode.vcf"
filepath <- paste0("/local/workdir/yc2644/CV_CT_array/vcf_n1601/", filename)
subset_pruned_vcf <- read.vcfR(filepath, verbose = FALSE)
subset_pruned_genind <- vcfR2genind(subset_pruned_vcf)

subset_pop_info <- pop_info <- read_tsv("../vcf_n1601/pop_grp_info_n1601.tsv", col_names = T) %>% 
  filter(Ind %in% readLines("../vcf_n1601/sample_id_n1601_contrast.txt"))

subset_pop_info$Pop = factor(subset_pop_info$Pop)
subset_pop_info$Location = factor(subset_pop_info$Location)
subset_pop_info$Grp = factor(subset_pop_info$Grp)
subset_pop_info$Type = factor(subset_pop_info$Type)

strata(subset_pruned_genind) <- subset_pop_info
subset_pruned_genind
ploidy(subset_pruned_genind) # all 2
tab(subset_pruned_genind)[1:5, 1:5] # confirm explicit allele dosages
# missing values trigger the warning; just set within=FALSE
amova.result1 <- poppr.amova(subset_pruned_genind, ~Location/Type, nperm=999,
                             within=FALSE)
amova.test1 <- randtest(amova.result1,nrepet = 999) # Test for significance

amova.result2 <- poppr.amova(subset_pruned_genind, ~Type/Location, nperm=999,
                             within=FALSE)
amova.test2 <- randtest(amova.result2,nrepet = 999) # Test for significance

# We expect variations within samples to give the greatest amount of variation for populations that are not significantly differentiated. 
# Sigma represents the variance, σ, for each hierarchical level and 
# to the right is the percent of the total.
# $statphi provides the ϕ, population differentiation statistics. 
# These are used to test hypotheses about population differentiation. 
# We would expect a higher ϕ
# statistic to represent a higher amount of differentiation.

# # subset Fst -------------------------------------------------------------
# subset_pruned_genind
# 
# subset_pop_info <- subset_pop_info %>% mutate(stage_yr=paste(stage, yr, sep = "_"))
# table(subset_pop_info$stage_yr)
# subset_pop_info$stage_yr = factor(subset_pop_info$stage_yr)
# 
# subset_pruned_genind@pop <- subset_pop_info$stage_yr
# subset_pruned_genind@pop <- subset_pop_info$yr
# subset_pruned_genind@pop <- subset_pop_info$stage
# 
# subset_pruned_genind
# # Weir and Cockerham's estimate
# wc(subset_pruned_genind)
# # $FST
# # [1] 0.0004552227
# # $FIS
# # [1] 0.1174404
# 
# # Pairwise Fst
# fst_pruned_subset <- genet.dist(subset_pruned_genind, method = "WC84") 

