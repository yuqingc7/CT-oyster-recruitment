setwd("/local/workdir/yc2644/CV_CT_array/vcf_n1601")

# see https://github.com/clairemerot/Intro2PopGenomics/blob/master/3.2.3/PCA/script_PCA_from_vcf.R for original code
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SeqArray")
library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/
library(tidyverse)
library(ggplot2)
library(SeqArray)
library(ggrepel)
library(RColorBrewer)

VCF="geno_n1601_fp_noinvers.vcf"

# metadata
meta <- read_tsv("pop_grp_info_n1601.tsv", col_names = T) %>% 
  mutate(
    Geo2 = case_when(
      Location %in% c("East",
                      "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "West CT",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East CT",
      TRUE ~ Grp
    )) 

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
table(tab_meta$Type)
table(tab_meta$Grp)
table(tab_meta$Geo2)

color = c("#8DA0CB", #AQ_1
          "#C77CFF", #AQ_2
          "#00BFC4", #AQ_3
          "#FF61C3", #AQ_4
          "#3BA4F4", #AQ_5
          "#3D2C7D", #AQ_6
          #"brown3", #CT
          #"#FC8D62", #ER_1
          "#911eb4","#66C2A5","#ffe764") #HR_4

p <- ggplot(tab_meta, aes(x = EV1, y = EV2,color=Geo2)) + 
  geom_point(size=3, aes(shape=Grp)) +
  scale_color_manual(values = color, name="Population",
                     #labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6","CT", "East", "Hudson")
                     ) +
  scale_fill_manual(values = color,name="Population",
                    #labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6","CT", "East", "Hudson")
                    ) +
  scale_shape_manual(values = c(16,16,16,16,16,16,10,10,10),name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "CT", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  coord_cartesian(xlim=c(-0.035,0.12))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  stat_ellipse(geom = "polygon",aes(fill = Grp),#type="norm",
               level = 0.95, alpha=0.05)+
  guides(color = "none",fill="none",shape="none")
p

tab_meta %>% 
  filter(Geo2 %in% c("East CT","West CT")) %>% 
  dplyr::select(-EV3,-EV4,-EV5,-EV6) %>% 
  filter(EV1>0.01)

tab_meta %>% 
  filter(Geo2 %in% c("East CT","West CT")) %>% 
  dplyr::select(-EV3,-EV4,-EV5,-EV6) %>% 
  filter(EV2>0.01)


ggplot(tab_meta, aes(x = EV3, y = EV4,color=Geo2)) + 
  geom_point(size=3, aes(shape=Grp)) +
  scale_color_manual(values = color, name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "CT", "East", "Hudson")) +
  scale_fill_manual(values = color,name="Population",
                    labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                               "CT", "East", "Hudson")) +
  scale_shape_manual(values = c(16,16,16,16,16,16,10,10,10),name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "CT", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC3 (",round(pc.percent[3],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC4 (",round(pc.percent[4],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.035,0.12))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  stat_ellipse(geom = "polygon",aes(fill = Grp),#type="norm",
               level = 0.95, alpha=0.05)+
 guides(color = "none",fill="none", shape="none")


# ref subset --------------------------------------------------------------

VCF="geno_HR_AQ_n284_fp_noinvers.vcf"

# metadata
meta <- read_tsv("pop_grp_info_n1601.tsv", col_names = T) %>% 
  filter(Ind %in% readLines("sample_id_n284_HR_AQ.txt"))

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
tab <- data.frame(Ind = readLines("sample_id_n284_HR_AQ.txt"),
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
table(tab_meta$Type)
table(tab_meta$Grp)

color = c("#8DA0CB", #AQ_1
          "#C77CFF", #AQ_2
          "#00BFC4", #AQ_3
          "#FF61C3", #AQ_4
          "#3BA4F4", #AQ_5
          "#3D2C7D", #AQ_6
          "#66C2A5") #HR_4

ggplot(tab_meta, aes(x = EV1, y = EV2,color=Grp)) + 
  geom_point(size=3, aes(shape=Grp)) +
  scale_color_manual(values = color, name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "Hudson")) +
  scale_fill_manual(values = color,name="Population",
                    labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                               "Hudson")) +
  scale_shape_manual(values = c(16,16,16,16,16,16,10,10,10),name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "Hudson")) +
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
        legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  stat_ellipse(geom = "polygon",aes(fill = Grp),#type="norm",
               level = 0.95, alpha=0.05)+
  guides(color = "none",fill="none",shape="none")

ggplot(tab_meta, aes(x = EV3, y = EV4,color=Grp)) + 
  geom_point(size=3, aes(shape=Grp)) +
  scale_color_manual(values = color, name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "Hudson")) +
  scale_fill_manual(values = color,name="Population",
                    labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                               "Hudson")) +
  scale_shape_manual(values = c(16,16,16,16,16,16,10,10,10),name="Population",
                     labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
                                "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC3 (",round(pc.percent[3],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC4 (",round(pc.percent[4],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.035,0.12))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  stat_ellipse(geom = "polygon",aes(fill = Grp),#type="norm",
               level = 0.95, alpha=0.05)+
  guides(color = "none",fill="none",shape="none")

# wild subset -------------------------------------------------------------

# wild subset
VCF="geno_n1367_fp_noinvers.recode.vcf"

meta <- read_tsv("pop_grp_info_n1601.tsv", col_names = T) %>% 
  filter(Ind %in% readLines("sample_id_n1601_wildsubset_n1367.txt")) %>% 
  mutate(
    Geo3 = case_when(
      Location %in% c("East",
                      "TOC","GOO","NRH","LOP","SAB") ~ "West",
      Location %in% c("HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "Central",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East",
      TRUE ~ Location
    )) %>%
  mutate(Geo3=factor(Geo3, levels=c("Hudson","Western", "Central", "Eastern"))) %>% 
  mutate(
    Geo2 = case_when(
      Location %in% c("East",
                      "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "West LIS",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East LIS",
      TRUE ~ Location
    )) %>% 
  mutate(Geo2=factor(Geo2, levels=c("Hudson","West LIS", "East LIS")))

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
table(tab_meta$Type)
table(tab_meta$Geo2)

x_zoom <- c(-0.045, 0.045)
y_zoom <- c(-0.05, 0.10)
ggplot(tab_meta, aes(x = EV1, y = EV2, color=Type)) + 
  geom_point(size=3, aes(shape=Type), alpha=0.8, stroke=0.5) +
  scale_color_manual(values = c("brown3","brown3","#FC8D62","#66C2A5"), name="Type",
                     labels = c("CT Adult", "CT Spat", "ER", "HR")) +
  # scale_fill_manual(values = color,name="Population",
  #                   labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                              "CT", "East", "Hudson")) +
  scale_shape_manual(values = c(1,2,17,17), name="Type",
                     labels = c("CT Adult", "CT Spat", "ER", "HR")) +
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
        legend.position = c(0.8, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  annotate("rect",
           xmin = x_zoom[1], xmax = x_zoom[2],
           ymin = y_zoom[1], ymax = y_zoom[2],
           color = "black", fill = NA, linewidth = 0.7, linetype = "dashed")

  # stat_ellipse(geom = "polygon",aes(fill = Type),#type="norm",
  #              level = 0.95, alpha=0.05)+
  # guides(color = "none",fill="none",shape="none")


# "#da73e6","lightgreen","#ffe764"
ggplot(tab_meta, aes(x = EV1, y = EV2, color=Geo2)) + 
  geom_point(size=3, aes(shape=Geo2)) +
  scale_color_manual(values = c("#66C2A5","#ffe764","#911eb4"),name="",
                     labels=c("HR", "West LIS", "East LIS")) +
  # scale_fill_manual(values = color,name="Population",
  #                   labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                              "CT", "East", "Hudson")) +
  scale_shape_manual(values = c(10,10,10),name="",labels=c("HR", "West LIS", "East LIS")) +
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
        panel.background = element_rect(colour = NA, size=0.5),
        legend.position = c(0.8, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()) +
  #labs(color = "") +
  annotate("rect",
           xmin = x_zoom[1], xmax = x_zoom[2],
           ymin = y_zoom[1], ymax = y_zoom[2],
           color = "black", fill = NA, linewidth = 0.7, linetype = "dashed")

# zoomed in
ggplot(tab_meta, aes(x = EV1, y = EV2, color=Type)) + 
  geom_point(size=3, aes(shape=Type), alpha=0.8, stroke=0.5) +
  scale_color_manual(values = c("brown3","brown3","#FC8D62","#66C2A5"), name="Type",
                     labels = c("CT Adult", "CT Spat", "ER", "HR")) +
  # scale_fill_manual(values = color,name="Population",
  #                   labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                              "CT", "East", "Hudson")) +
  scale_shape_manual(values = c(1,2,17,17), name="Type",
                     labels = c("CT Adult", "CT Spat", "ER", "HR")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  coord_cartesian(xlim=c(x_zoom[1],x_zoom[2]),
                  ylim=c(y_zoom[1],y_zoom[2]))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        legend.position = c(0.75, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  guides(color = "none",fill="none",shape="none")

# "#da73e6","lightgreen","#ffe764"
ggplot(tab_meta, aes(x = EV1, y = EV2, color=Geo2)) + 
  geom_point(size=3, aes(shape=Geo2)) +
  scale_color_manual(values = c("#66C2A5","#ffe764","#911eb4"),name="",
                     labels=c("HR", "West CT", "East CT")) +
  # scale_fill_manual(values = color,name="Population",
  #                   labels = c("AQ_1", "AQ_2", "AQ_3", "AQ_4","AQ_5", "AQ_6",
  #                              "CT", "East", "Hudson")) +
  scale_shape_manual(values = c(10,10,10),name="", labels=c("HR", "West CT", "East CT")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  coord_cartesian(xlim=c(x_zoom[1],x_zoom[2]),
                  ylim=c(y_zoom[1],y_zoom[2]))+
  #coord_cartesian(xlim=c(-0.035,0.12))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = NA, size=0.5),
        legend.position = c(0.8, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()) +
  guides(color = "none",fill="none",shape="none")
  #labs(color = "") 

# explore -----------------------------------------------------------------

library(ggforce)
library(colorspace)
p <- ggplot(tab_meta, aes(x = EV1, y = EV2,color=Location)) + 
  geom_point(size=3, aes(shape=Type)) +
  scale_color_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_fill_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_shape_manual(values = c(16,10,10,10), name="Type",
                     labels = c("CT Adult", "CT Spat", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.1,0.05))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        #legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  # stat_ellipse(geom = "polygon",aes(fill = Location),#type="norm",
  #              level = 0.95, alpha=0.05, label= "Location") +
  #geom_mark_ellipse(aes(label = Pop), alpha=0.01) +
  #guides(color = "none",fill="none",shape="none")
  geom_text(data = subset(tab_meta, grepl("YBD", Ind)), 
            aes(label = Ind), hjust = 1.2, vjust = 1.2, size = 3) 
# geom_text(data = subset(tab_meta, EV1 < -0.025), 
#         aes(label = Ind), hjust = 1.2, vjust = 1.2, size = 3)

p

ggplot(data.frame(1:32,pc.percent[1:32]), aes(x = 1:32, y = pc.percent[1:32])) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Explained Variance") +
  scale_x_continuous(breaks = seq(1, length(pc.percent), 1)) +  
  theme_minimal()

ggplot(tab_meta, aes(x = EV3, y = EV4,color=Location)) + 
  geom_point(size=3, aes(shape=Type)) +
  scale_color_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_fill_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_shape_manual(values = c(16,10,10,10), name="Type",
                     labels = c("CT Adult", "CT Spat", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC3 (",round(pc.percent[3],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC4 (",round(pc.percent[4],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.1,0.05))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        #legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  # stat_ellipse(geom = "polygon",aes(fill = Location),#type="norm",
  #              level = 0.95, alpha=0.05, label= "Location") +
  #geom_mark_ellipse(aes(label = Pop), alpha=0.01) +
  #guides(color = "none",fill="none",shape="none")
  geom_text(data = subset(tab_meta, grepl("YBD", Ind)), 
            aes(label = Ind), hjust = 1.2, vjust = 1.2, size = 3) 

ggplot(tab_meta, aes(x = EV5, y = EV6,color=Location)) + 
  geom_point(size=3, aes(shape=Type)) +
  scale_color_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_fill_manual(values = qualitative_hcl(23, palette = "Dark 3"), name="Location") +
  scale_shape_manual(values = c(16,10,10,10), name="Type",
                     labels = c("CT Adult", "CT Spat", "East", "Hudson")) +
  # scale_shape_manual(values=c(16,10), name="",
  #                    labels = c("Aquaculture", "Native")) +
  scale_x_continuous(paste("PC5 (",round(pc.percent[5],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC6 (",round(pc.percent[6],3),"%",")",sep=""))+ 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  theme_bw()+
  #coord_cartesian(xlim=c(-0.1,0.05))+
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        #legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0)) +
  labs(color = "") +
  # stat_ellipse(geom = "polygon",aes(fill = Location),#type="norm",
  #              level = 0.95, alpha=0.05, label= "Location") +
  #geom_mark_ellipse(aes(label = Pop), alpha=0.01) +
  #guides(color = "none",fill="none",shape="none")
  geom_text(data = subset(tab_meta, grepl("YBD", Ind)), 
            aes(label = Ind), hjust = 1.2, vjust = 1.2, size = 3) 
