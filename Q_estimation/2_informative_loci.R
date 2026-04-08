setwd("/local/workdir/yc2644/CV_CT_array/vcf_n1601")

library(hierfstat)
library(vcfR)
library(adegenet)

vcf <- read.vcfR("geno_HR_AQ_n284_fp_noinvers.vcf")
genind <- vcfR2genind(vcf)

# Fst 2 pops --------------------------------------------------------------------

fst_values <- read_tsv("HR_n50_AQ_n234_fst.weir.fst")
# CHROM   POS     WEIR_AND_COCKERHAM_FST
# 1       30115   -0.00438485
# 1       38807   0.0516085
# 1       44310   0.212624
fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")

top1 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.99)
top5 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.95)
top10 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.90)
top15 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.85)
top20 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.80)

top_fst_1 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top1] #386
top_fst_5 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top5] #1926
top_fst_10 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top10] #3852
top_fst_15 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top15] #5777
top_fst_20 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top20] #7705


# DAPC (2 pops) --------------------------------------------------------------------

# Read popmap: assuming 2 columns, sample ID and population``
popmap <- read_csv("../lists/HR_AQ_n284_map.txt", col_names = T)

# Match samples to genind object
genind$pop <- factor(popmap$pop_id[match(indNames(genind), popmap$sample_id)])

# tab() converts genind to numeric matrix
geno_matrix <- tab(genind, NA.method = "mean")  # replace NAs by mean allele frequency per SNP

# Run PCA (retain enough PCs to capture most variation)
#pca_res <- dudi.pca(geno_matrix, scannf = FALSE, nf = 50)  # nf = number of PCs

# Run DAPC
#dapc_test <- dapc(geno_matrix, genind$pop, n.pca = 50, 
#                  n.da = length(unique(genind$pop)) - 1)
optim_res <- optim.a.score(dapc_test)  # suggests optimal number of PCs

dapc_res <- dapc(geno_matrix, genind$pop, n.pca = 7, 
                 n.da = length(unique(genind$pop)) - 1)

# Plot
scatter(dapc_res, scree.da = F, posi.da="bottomright", cstar=1, clabel=0.75)
scatter(
  dapc_res,
  #col = colors,           # color by population
  scree.da = FALSE,       # do not plot scree of DFs
  posi.da = "bottomright",# label position
  cstar = 0,              # no lines from individuals to centroid
  clabel = 0.75           # scaling of population labels
)
loadingplot(dapc_res$var.contr, threshold=0.01, lab.jitter = 1.2)

# Contribution per SNP
snp_contrib <- dapc_res$var.contr

# If multiple DAs, sum contributions across all DFs
snp_total_contrib <- rowSums(snp_contrib)

# Rank SNPs by contribution
snp_ranked <- names(sort(snp_total_contrib, decreasing = TRUE))
snp_ranked_clean <- sub("\\..*$", "", snp_ranked)

top_dapc_1 <- unique(snp_ranked_clean[1:round(0.01 * length(snp_ranked_clean))]) #385
top_dapc_5 <- unique(snp_ranked_clean[1:round(0.05 * length(snp_ranked_clean))]) #1926
top_dapc_10 <- unique(snp_ranked_clean[1:round(0.1 * length(snp_ranked_clean))]) #3852
top_dapc_15 <- unique(snp_ranked_clean[1:round(0.15 * length(snp_ranked_clean))]) #5777
top_dapc_20 <- unique(snp_ranked_clean[1:round(0.20 * length(snp_ranked_clean))]) #7703


# Top SNPs meeting ALL criteria -------------------------------------------

selected_1 <- intersect(top_fst_1, top_dapc_1) #103
selected_5 <- intersect(top_fst_5, top_dapc_5) #905
selected_10 <- intersect(top_fst_10, top_dapc_10) #2186
selected_15 <- intersect(top_fst_15, top_dapc_15) #3592
selected_20 <- intersect(top_fst_20, top_dapc_20) #5089

writeLines(selected_5, "../SNP_selected_5.list")



# # Fst 7 pops --------------------------------------------------------------------

# fst_files <- list.files(pattern = "HR_n50_.*_fst\\.weir\\.fst$")[-1]
# fst_files

# # function to read fst file and extract top 25% SNPs
# get_top25_snps <- function(file) {
#   fst_values <- read_tsv(file,col_names = TRUE) %>% 
#     drop_na(WEIR_AND_COCKERHAM_FST)
#   fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")
  
#   top25 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.75, na.rm = TRUE)
#   top_fst_25 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top25]
  
#   return(top_fst_25)
# }

# # function to read fst file and extract top 30% SNPs
# get_top30_snps <- function(file) {
#   fst_values <- read_tsv(file,col_names = TRUE) %>% 
#     drop_na(WEIR_AND_COCKERHAM_FST)
#   fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")
  
#   top30 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.7, na.rm = TRUE)
#   top_fst_30 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top30]
  
#   return(top_fst_30)
# }

# get_top10_snps <- function(file) {
#   fst_values <- read_tsv(file,col_names = TRUE) %>% 
#     drop_na(WEIR_AND_COCKERHAM_FST)
#   fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")
  
#   top10 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.90, na.rm = TRUE)
#   top_fst_10 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top10]
  
#   return(top_fst_10)
# }

# get_top5_snps <- function(file) {
#   fst_values <- read_tsv(file,col_names = TRUE) %>% 
#     drop_na(WEIR_AND_COCKERHAM_FST)
#   fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")
  
#   top5 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.95, na.rm = TRUE)
#   top_fst_5 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top5]
  
#   return(top_fst_5)
# }

# get_top1_snps <- function(file) {
#   fst_values <- read_tsv(file,col_names = TRUE) %>% 
#     drop_na(WEIR_AND_COCKERHAM_FST)
#   fst_values$SNP_ID <- paste(fst_values$CHROM, fst_values$POS, sep = "_")
  
#   top1 <- quantile(fst_values$WEIR_AND_COCKERHAM_FST, 0.99, na.rm = TRUE)
#   top_fst_1 <- fst_values$SNP_ID[fst_values$WEIR_AND_COCKERHAM_FST >= top1]
  
#   return(top_fst_1)
# }

# # apply to all files
# top25_list <- lapply(fst_files, get_top25_snps)
# names(top25_list) <- fst_files
# top30_list <- lapply(fst_files, get_top30_snps)
# names(top30_list) <- fst_files
# top10_list <- lapply(fst_files, get_top10_snps)
# names(top10_list) <- fst_files
# top5_list <- lapply(fst_files, get_top5_snps)
# names(top5_list) <- fst_files
# top1_list <- lapply(fst_files, get_top1_snps)
# names(top1_list) <- fst_files

# # find intersection across all datasets
# intersect_fst_top25 <- Reduce(intersect, top25_list)
# length(intersect_fst_top25) #124
# intersect_fst_top30 <- Reduce(intersect, top30_list)
# length(intersect_fst_top30) #206

# union_fst_top10 <- Reduce(union, top10_list)
# length(union_fst_top10)
# union_fst_top5 <- Reduce(union, top5_list)
# length(union_fst_top5)
# union_fst_top1 <- Reduce(union, top1_list)
# length(union_fst_top1)

# # DAPC (7 pops) --------------------------------------------------------------------

# # Read popmap: assuming 2 columns, sample ID and population``
# popmap <- read_csv("../lists/HR_AQ_n284_map2.txt", col_names = T)

# # Match samples to genind object
# genind$pop <- factor(popmap$pop_id[match(indNames(genind), popmap$sample_id)])

# # tab() converts genind to numeric matrix
# geno_matrix <- tab(genind, NA.method = "mean")  # replace NAs by mean allele frequency per SNP

# # Run PCA (retain enough PCs to capture most variation)
# #pca_res <- dudi.pca(geno_matrix, scannf = FALSE, nf = 50)  # nf = number of PCs

# # Run DAPC
# dapc_test <- dapc(geno_matrix, genind$pop, n.pca = 100, 
#                   n.da = length(unique(genind$pop)) - 1)
# optim_res <- optim.a.score(dapc_test)  # suggests optimal number of PCs

# dapc_res <- dapc(geno_matrix, genind$pop, n.pca = 11, 
#                  n.da = length(unique(genind$pop)) - 1)

# color = c("#8DA0CB", #AQ_1
#           "#C77CFF", #AQ_2
#           "#00BFC4", #AQ_3
#           "#FF61C3", #AQ_4
#           "#3BA4F4", #AQ_5
#           "#3D2C7D", #AQ_6
#           "#66C2A5") #HR_4

# # Plot
# df_percent <- 100 * dapc_res$eig / sum(dapc_res$eig)

# # For DF1 and DF2 (if available)
# df1_perc <- round(df_percent[1], 1)
# df2_perc <- if(length(df_percent) > 1) round(df_percent[2], 1) else NA

# scatter(dapc_res,cex = 2, legend = F,col=color,cstar=1,
#         clabel = F, posi.leg = "right", scree.pca = F,posi.da="bottom",
#         cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

# loadingplot(dapc_res$var.contr, threshold=0.01, lab.jitter = 1.2)

# # Contribution per SNP
# snp_contrib <- dapc_res$var.contr

# # If multiple DAs, sum contributions across all DFs
# snp_total_contrib <- rowSums(snp_contrib)

# # Rank SNPs by contribution
# snp_ranked <- names(sort(snp_total_contrib, decreasing = TRUE))
# snp_ranked_clean <- sub("\\..*$", "", snp_ranked)

# top_dapc7_1 <- unique(snp_ranked_clean[1:round(0.01 * length(snp_ranked_clean))])
# top_dapc7_5 <- unique(snp_ranked_clean[1:round(0.05 * length(snp_ranked_clean))])
# top_dapc7_10 <- unique(snp_ranked_clean[1:round(0.1 * length(snp_ranked_clean))])
# top_dapc7_15 <- unique(snp_ranked_clean[1:round(0.15 * length(snp_ranked_clean))])
# top_dapc7_20 <- unique(snp_ranked_clean[1:round(0.20 * length(snp_ranked_clean))])
# top_dapc7_25 <- unique(snp_ranked_clean[1:round(0.25 * length(snp_ranked_clean))])
# top_dapc7_30 <- unique(snp_ranked_clean[1:round(0.30 * length(snp_ranked_clean))])

