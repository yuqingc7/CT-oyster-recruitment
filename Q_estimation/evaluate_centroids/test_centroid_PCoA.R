
setwd("/local/workdir/yc2644/CV_CT_array/centroid_sim/check_centroids")

library(ape)
library(ggplot2)

# FST matrix (4 populations: AQ, HR, in_silico_AQ, in_silico_wild)
fst <- matrix(c(
  0,          0.029925, -0.0013262, 0.034924,
  0.029925,   0,         0.031237,  -0.0052743,
  -0.0013262,  0.031237,  0,         0.03613,
  0.034924,  -0.0052743, 0.03613,   0
), nrow=4, byrow=TRUE)

rownames(fst) <- colnames(fst) <- c("AQ","HR","in_silico_AQ","in_silico_wild")

# Make sure it's symmetric
fst <- (fst + t(fst))/2
diag(fst) <- 0
fst

pcoa_res <- pcoa(fst)

# Extract coordinates for first two axes
pcoa_df <- data.frame(
  Population = rownames(fst),
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2],
  Group = c("Empirical","Empirical","In_Silico","In_Silico")
)

ggplot(pcoa_df, aes(x=Axis1, y=Axis2, label=Population, color=Group)) +
  geom_point(size=4) +
  geom_text(vjust=-0.8) +
  scale_color_manual(values=c("Empirical"="#1f78b4","In_Silico"="#e31a1c")) +
  theme_minimal() +
  labs(
    x=paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
    y=paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)"),
    title="PCoA of pairwise FST"
  ) +
coord_cartesian(ylim = c(-0.01, 0.01),xlim = c(-0.03, 0.03))


# 6 AQs -------------------------------------------------------------------

# Step 1: Read grep output
lines <- readLines("mean_FST_all.txt")

# Step 2: Extract population names and FST
pop1 <- sapply(strsplit(lines, "_vs_"), `[`, 1)
pop2 <- sapply(strsplit(sapply(strsplit(lines, "_vs_"), `[`, 2), ".log"), `[`, 1)
fst_vals <- as.numeric(sub(".*: ", "", lines))

# Step 3: Unique populations
pops <- unique(c(pop1, pop2))
clean_pops <- gsub("geno_|_samples", "", pops)

n <- length(pops)

# Step 4: Build empty matrix
fst_mat <- matrix(0, nrow=n, ncol=n)
rownames(fst_mat) <- colnames(fst_mat) <- pops

# Step 5: Fill symmetric matrix
for(i in seq_along(fst_vals)){
  fst_mat[pop1[i], pop2[i]] <- fst_vals[i]
  fst_mat[pop2[i], pop1[i]] <- fst_vals[i]
}

# Optional: set negative FST to zero
#fst_mat[fst_mat < 0] <- 0

# Step 6: Run PCoA
pcoa_res <- pcoa(fst_mat)

# Step 7: Prepare plotting dataframe
pcoa_df <- data.frame(
  Population = clean_pops, # need to match rownames(fst_mat)
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2]
)

# Step 8: Assign groups for coloring
group <- ifelse(grepl("in_silico", pcoa_df$Population), "In Silico",
                ifelse(grepl("HR", pcoa_df$Population), "HR Empirical", "AQ Empirical"))
pcoa_df$Group <- group

# Step 9: Pretty 2D PCoA plot
ggplot(pcoa_df, aes(x=Axis1, y=Axis2, color=Group, fill=Group)) +
  geom_point(aes(shape=Group),size=5, stroke=1) +
  geom_text(aes(label=Population), 
            size=4, 
            nudge_y = 0.02,   # move labels slightly up
            #nudge_x = 0.02,   # move labels slightly right
            check_overlap = T) +  # removes labels that still overlap
  stat_ellipse(alpha=0.2, geom="polygon", color=NA) +
  scale_shape_manual(values=c("AQ Empirical"=18,"HR Empirical"=18,
                              "In Silico"=8)) +
  scale_color_manual(values=c("AQ Empirical"="#1f78b4","HR Empirical"="#66C2A5",
                              "In Silico"="#e31a1c")) +
  scale_fill_manual(values=c("AQ Empirical"="#1f78b4","HR Empirical"="#66C2A5",
                             "In Silico"="#e31a1c")) +
  theme_minimal(base_size=14) +
  theme(legend.position="") +
  labs(
    x=paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1]*100,1), "%)"),
    y=paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2]*100,1), "%)"),
    title="PCoA of pairwise FST (Selected 905 SNPs)"
  ) +
  coord_cartesian(ylim = c(-0.25, 0.25), xlim = c(-0.25, 0.25))




# try PCA -----------------------------------------------------------------

library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/
library(tidyverse)
library(ggplot2)
library(SeqArray)
library(ggrepel)
library(RColorBrewer)

VCF="geno_allrefs_wild_AQ_fp_noinvers.recode.vcf"

# metadata
meta <- read_tsv("geno_allrefs_wild_AQ_samples.txt", col_names = F) %>% 
  mutate(Pop = sub("-.*", "", X1)) %>% 
  rename(Ind=X1)

# vcf
vcf.fn <- paste0("./",VCF)
# VCF => GDS
showfile.gds(closeall=TRUE)
snpgdsVCF2GDS(vcf.fn, "pruned.noinvers.gds", method="biallelic.only")
# showfile.gds(closeall=TRUE)
# summary
snpgdsSummary("pruned.noinvers.gds")
# Open the GDS file
genofile <- snpgdsOpen("pruned.noinvers.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(Ind = meta$Ind,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
dim(tab)
dim(meta)

tab_meta <- right_join(meta,tab, by="Ind")
table(tab_meta$Pop)

# color = c("#8DA0CB", #AQ_1
#           "#C77CFF", #AQ_2
#           "#00BFC4", #AQ_3
#           "#FF61C3", #AQ_4
#           "#3BA4F4", #AQ_5
#           "#3D2C7D", #AQ_6
#           "brown3", #CT
#           "#FC8D62", #ER_1
#           "#66C2A5") #HR_4

ggplot(tab_meta, aes(x = EV1, y = EV2,color=Pop)) + 
  geom_point(size=3) +
  # scale_color_manual(values = color, name="Population",
  #                    labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                               "CT", "East", "Hudson")) +
  # scale_fill_manual(values = color,name="Population",
  #                   labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                              "CT", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.035,0.12))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        legend.position = "top",  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") 
  # stat_ellipse(geom = "polygon",aes(fill = Pop),#type="norm",
  #              level = 0.95, alpha=0.05)+
  #guides(color = "none",fill="none",shape="none")
