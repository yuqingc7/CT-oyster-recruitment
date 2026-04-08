setwd("/local/workdir/yc2644/CV_CT_array/structure")

library(RColorBrewer)
library(tidyverse)
library(pophelper)
library(viridis)
library(gridExtra)

### Read in data------------------------------------------------------------
# Q tables
n1601_k5 <- readQ("./n1601_cor/str_K5_rep1_f", filetype="auto")$str_K5_rep1_f # best K
n1601_k4 <- readQ("./n1601_cor/str_K4_rep1_f", filetype="auto")$str_K4_rep1_f
n1601_k6 <- readQ("./n1601_cor/str_K6_rep1_f", filetype="auto")$str_K6_rep1_f

#n1601_k5 <- readQ("./n1601_snp225_cor/str_K5_rep1_f", filetype="auto")$str_K5_rep1_f # best K

n1601 <- n1601_k5
# ADMIXTURE result
#tbl <- read.table("/workdir/yc2644/CV_CT_array/ADMIXTURE/10kSNP.geno_n1601_fp_noinvers.6.Q")
# n1601 <- tbl %>% 
#   rename(Cluster1=V1,Cluster2=V2,Cluster3=V3,Cluster4=V4,
#          Cluster5=V5,Cluster6=V6)

# metadata
grps <- read_tsv("../vcf_n1601/pop_grp_info_n1601.tsv") %>% 
  mutate(
    Geo2 = case_when(
      Location %in% c("East",
                      "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "West CT",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East CT",
      TRUE ~ Grp
    )) 

dim(grps)

# combine sample IDs to qlist
str_n1601 <- cbind(grps, n1601) 

#write_tsv(str_n1601,"str_n1601_k5.txt", col_names = T)

### Summarize data------------------------------------------------------------

str_n1601 %>% 
  filter(Grp %in% c("CT","East","Hudson")) %>% 
ggplot(aes(x = 1-Cluster4, fill = Grp)) +
 geom_histogram(#aes(y = after_stat(density)),
                 binwidth = 0.005,
                 alpha = 1,
                 position = "identity") +
  # geom_density(color = "black", 
  #                alpha = 0.6, 
  #                position = "identity") +
  geom_vline(xintercept=zero_HR, color="blue", linetype="dashed") +
  facet_wrap(~Grp, scales = "free_y",ncol=1) +
  scale_fill_manual(values=c("brown3", #CT
                             "#FC8D62", #ER_1
                             "#66C2A5")) +
  theme_bw()+ 
  #labs(x = "Admixture Level (Q)", y = "Density")+
  labs(x = "Admixture Level (Q)", y = "Count") +
  theme(axis.text=element_text(size=12),
        #text = element_text(size=18,family="Times"),
        text = element_text(size=15),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5),
        #legend.position = c(0.15, 0.8),  # Adjust these values to move the legend
        legend.margin = margin(0, 0, 0, 0))

mean(1-subset(str_n1601, Grp %in% c("AQ_1","AQ_2","AQ_3","AQ_5","AQ_6"))$Cluster4)
mean(1-subset(str_n1601, Grp %in% c("Hudson"))$Cluster4)

table(str_n1601$Grp)
zero_HR <- quantile(1-subset(str_n1601, Grp=="Hudson")$Cluster4, probs=0.90) #0.008
#quantile(1-subset(str_n1601, Grp=="Hudson")$Cluster4, probs=0.95)

HR <- str_n1601 %>% filter(Grp %in% c("Hudson")) #50
unadmixed_HR <- HR$Ind[HR$Cluster4 >= 1-zero_HR]
length(unadmixed_HR) #46
#writeLines(unadmixed_HR, "../lists/unadmixed_new_HR_n46.txt")

ggplot(subset(str_n1601, Grp=="Hudson"), aes(x = 1-Cluster4)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  geom_vline(xintercept = zero_HR, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Admixture Level (Q)", y = "Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

# proportion of pop being admixed or not d
dim(str_n1601 %>% filter(Grp %in% c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6"))) #234

CT_ER <- str_n1601 %>% filter(Grp %in% c("CT", "East")) #1317
#writeLines(CT_ER$Ind, "../lists/CT_ER_new_n1317.txt")

admixed <- CT_ER$Ind[CT_ER$Cluster4 < 1-zero_HR]
length(admixed) #722
722/1317
#writeLines(admixed, "../lists/admixed_new_CT_ER_n722.txt")

admixed_CT_spat <- CT_ER %>% filter(Cluster4 < 1-zero_HR, Type=="CT_spat") #467
admixed_CT_adult <- CT_ER %>% filter(Cluster4 < 1-zero_HR, Type=="CT_adult") #228
# writeLines(admixed_CT_spat$Ind, "../lists/admixed_new_CTspat_n467.txt")
# writeLines(admixed_CT_adult$Ind, "../lists/admixed_new_CTadult_n228.txt")

# exclude >F1
admixed_CT_spat_0.25 <- admixed_CT_spat %>% filter(Cluster4 > 0.750) #460
admixed_CT_adult_0.25 <- admixed_CT_adult %>% filter(Cluster4 > 0.750) #224
# writeLines(admixed_CT_spat_0.25$Ind, "../lists/admixed0.25_new_CTspat_n460.txt")
# writeLines(admixed_CT_adult_0.25$Ind, "../lists/admixed0.25_new_CTadult_n224.txt")
admixed_CT_ER_0.25 <- CT_ER %>% filter(Cluster4 > 0.750) %>% filter(Cluster4 < 1-zero_HR) #710
# writeLines(admixed_CT_ER_0.25$Ind, "../lists/admixed0.25_new_CT_ER_n710.txt")

unadmixed_CT_ER <- CT_ER$Ind[CT_ER$Cluster4 >= 1-zero_HR]
length(unadmixed_CT_ER) #595
#writeLines(unadmixed_CT_ER, "../lists/unadmixed_new_CT_ER_n595.txt")

CT <- str_n1601 %>% filter(Grp %in% c("CT")) #1259
admixed_CT <- CT$Ind[CT$Cluster4 < 1-zero_HR]
length(admixed_CT) #695
695/1259 
#writeLines(admixed_CT, "../lists/admixed_new_CT_n695.txt")

unadmixed_CT <- CT$Ind[CT$Cluster4 >= 1-zero_HR]
length(unadmixed_CT) #564
564/1259
#writeLines(unadmixed_CT, "../lists/unadmixed_new_CT_n564.txt")

CT_noYBD <- str_n1601 %>% filter(Grp %in% c("CT")) %>% 
  filter(Pop != "YBD") #1225
admixed_CT_noYBD <- CT_noYBD$Ind[CT_noYBD$Cluster4 < 1-zero_HR]
length(admixed_CT_noYBD) #670
670/1225

ER <- str_n1601 %>% filter(Grp %in% c("East")) #58
admixed_ER <- ER$Ind[ER$Cluster4 < 1-zero_HR]
length(admixed_ER) #27
27/58

unadmixed_ER <- ER$Ind[ER$Cluster4 >= 1-zero_HR]
length(unadmixed_ER)
#writeLines(unadmixed_ER, "../lists/unadmixed_new_ER_n31.txt")


# mean admixture excluding unadmixed inds
1-mean((subset(CT_ER, Cluster4 < 1-zero_HR))$Cluster4) #0.05803186
1-mean((subset(ER, Cluster4 < 1-zero_HR))$Cluster4) #0.06803704
1-mean((subset(CT, Cluster4 < 1-zero_HR))$Cluster4) #0.05764317
1-mean((subset(CT_noYBD, Cluster4 < 1-zero_HR))$Cluster4) #0.05480597


## adult vs spat | site -----------------------------------------------------------

# focus on where adult-spat pairs are available
table(CT$Type)
table(CT$Location)
table(CT$Pop)
CT_pair <- CT %>% filter(Location %in% c("ASC","EHF","FEC","GOO",
                              "HDC","JAC","SPT","TOC")) #761
admixed_CT_pair <- CT_pair$Ind[CT_pair$Cluster4 < 1-zero_HR]
length(admixed_CT_pair) #396
396/761
1-mean((subset(CT_pair, Cluster4 < 1-zero_HR))$Cluster4) #0.05388636

table(CT_pair$Type)
CT_pair_s <- CT_pair %>% filter(Type %in% c("CT_spat")) #466
admixed_CT_pair_s <- CT_pair_s$Ind[CT_pair_s$Cluster4 < 1-zero_HR]
length(admixed_CT_pair_s) #250
250/466 #0.5364807
1-mean((subset(CT_pair_s, Cluster4 < 1-zero_HR))$Cluster4) #0.051236

CT_pair_a <- CT_pair %>% filter(Type %in% c("CT_adult")) #295
admixed_CT_pair_a <- CT_pair_a$Ind[CT_pair_a$Cluster4 < 1-zero_HR]
length(admixed_CT_pair_a) #146
146/295 #0.4949153
1-mean((subset(CT_pair_a, Cluster4 < 1-zero_HR))$Cluster4) #0.05842466


### Prepare for plotting -----------------------------------------------------
# set population order
table(str_n1601$Pop)
table(str_n1601$Location)
table(str_n1601$Grp)
table(str_n1601$Type)
table(str_n1601$Geo2)
pop_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6",  
               "Hudson","East",
               "TOC0423","TOC23a","GOO0423","GOO23a","NRH0423","LOP0423","SAB0423","HOR",
               "ASC0423","ASC23a", "LIW1",
               "SPT22s", "STP23a", "MIH0423","MCK0423","QUR23a", "EHF0523","EHF23a",
               "Cv5785","JAC0523","JAC23a","HDC0523","HDC23a",
               "FEC0423","FEC23a","HAM0423","BAK0423","LIW2","YBD")
loc_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6",  
               "Hudson","East",
               "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
               "SPT","MIH","MCK","QUR", "EHF","Cv5785","JAC","HDC",
               "FEC","HAM","BAK","LIW2","YBD")
grp_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6", "Hudson", "East", "CT")
str_n1601 <- str_n1601 %>%
  mutate(Location=factor(Location, levels = loc_order),
         Grp=factor(Grp, levels = grp_order),
         Pop=factor(Pop, levels = pop_order),
         Geo2=factor(Geo2, 
                     levels= c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6","Hudson","West CT", "East CT"))) %>%
  arrange(Geo2,Pop)  

### Final Plots  -----------------------------------------------------

# get vline intercepts
table(str_n1601$Geo2)
counts <- as.vector(table(str_n1601$Geo2))
# cumulative sum
cumulative <- cumsum(counts)
cumulative <- cumulative[1:length(counts)-1]

# ordered by geographic locations
str_n1601_geo <- str_n1601 %>% 
  mutate(id = 1:1601) %>% 
  pivot_longer(cols = c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5"),
               names_to = "Source", values_to = "Q")
str_n1601_geo

ggplot(str_n1601_geo, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title="K=5")+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

# ordered by admixture levels within each pop
table(str_n1601$Pop)
# get vline intercepts
table(str_n1601$Location)
new_counts <- as.vector(table(str_n1601$Location))
# cumulative sum
new_cumulative <- cumsum(new_counts)
new_cumulative
str_n1601_anspop <- str_n1601 %>% 
  arrange(Pop, desc(1-Cluster4)) %>% # update accordingly
  mutate(id = 1:1601) %>% 
  pivot_longer(cols = c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5"), 
               names_to = "Source", values_to = "Q") 
str_n1601_anspop

ggplot(str_n1601_anspop, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = new_cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title="")+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")
table(str_n1601$Pop)

# ordered by admixture levels within each location
str_n1601_ans <- str_n1601 %>% 
  arrange(Location, desc(1-Cluster4)) %>% # update accordingly
  mutate(id = 1:1601) %>% 
  pivot_longer(cols = c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5"), 
               names_to = "Source", values_to = "Q") 
str_n1601_ans

ggplot(str_n1601_ans, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title="")+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

# ordered by admixture levels within Grp
str_n1601_ansgrp <- str_n1601 %>% 
  arrange(Grp, desc(1-Cluster4)) %>% 
  mutate(id = 1:1601) %>% 
  pivot_longer(cols = c("Cluster1", "Cluster2","Cluster3","Cluster4","Cluster5"), 
               names_to = "Source", values_to = "Q") 
str_n1601_ansgrp

ggplot(str_n1601_ansgrp, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title="K=5")+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")


# Plot K=4-6 --------------------------------------------------------------
library(patchwork)

# metadata
grps <- read_tsv("../vcf_n1601/pop_grp_info_n1601.tsv") %>% 
  mutate(
    Geo2 = case_when(
      Location %in% c("East",
                      "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "West CT",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East CT",
      TRUE ~ Grp
    )) 
dim(grps)
loc_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6",  
               "Hudson","East",
               "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
               "SPT","MIH","MCK","QUR", "EHF","Cv5785","JAC","HDC",
               "FEC","HAM","BAK","LIW2","YBD")
grp_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6", "Hudson", "East", "CT")

n1601_k5 <- readQ("./n1601_cor/str_K5_rep1_f", filetype="auto")$str_K5_rep1_f # best K
n1601_k4 <- readQ("./n1601_cor/str_K4_rep1_f", filetype="auto")$str_K4_rep1_f
n1601_k6 <- readQ("./n1601_cor/str_K6_rep1_f", filetype="auto")$str_K6_rep1_f
n1601_k3 <- readQ("./n1601_cor/str_K3_rep1_f", filetype="auto")$str_K3_rep1_f
n1601_k2 <- readQ("./n1601_cor/str_K2_rep1_f", filetype="auto")$str_K2_rep1_f

n1601_k2 <- readQ("./n1601_snp225_cor/str_K2_rep1_f", filetype="auto")$str_K2_rep1_f
dim(n1601_k2)

prep_qdata <- function(qdf, grps, k) {
  df <- cbind(grps, qdf) %>% # combine sample IDs to Q
    mutate(Grp = factor(Grp, levels = grp_order),
           Location = factor(Location, levels = loc_order),
           Geo2=factor(Geo2, 
                       levels= c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6","Hudson","West CT", "East CT"))) %>%
    arrange(Geo2, Location)
  return(df)
}

# Prepare data for each K
df_k2 <- prep_qdata(n1601_k2, grps, 2)
df_k3 <- prep_qdata(n1601_k3, grps, 3)
df_k4 <- prep_qdata(n1601_k4, grps, 4)
df_k5 <- prep_qdata(n1601_k5, grps, 5)
df_k6 <- prep_qdata(n1601_k6, grps, 6)

# get vline intercepts
table(df_k4$Geo2)
counts <- as.vector(table(df_k4$Geo2))
# cumulative sum
cumulative <- cumsum(counts)
cumulative <- cumulative[1:length(counts)-1]


# ordered by geographic locations
df_k2 <- df_k2 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 2))
df_k3 <- df_k3 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 3))
df_k4 <- df_k4 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 4))
df_k5 <- df_k5 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 5))
df_k6 <- df_k6 %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 6))

ggplot(df_k2, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k2$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

ggplot(df_k3, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k3$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

ggplot(df_k4, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k4$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

ggplot(df_k6, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k6$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

ggplot(df_k5, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k5$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none") +
  coord_cartesian(c(1550,1601))

df_k2 <- df_k2 %>%
  # ordered by admixture levels within each pop
  arrange(Location, desc(1-Cluster1)) %>% # update accordingly
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 2))

df_k3 <- df_k3 %>%
  # ordered by admixture levels within each pop
  arrange(Location, desc(1-Cluster3)) %>% # update accordingly
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 3))

df_k4 <- df_k4 %>%
  # ordered by admixture levels within each pop
  arrange(Location, desc(1-Cluster3)) %>% # update accordingly
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 4))

df_k5 <- df_k5 %>%
  # ordered by admixture levels within each pop
  arrange(Location, desc(1-Cluster4)) %>% # update accordingly
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 5))

df_k6 <- df_k6 %>%
  # ordered by admixture levels within each pop
  arrange(Location, desc(1-Cluster2)) %>% # update accordingly
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "Source", values_to = "Q") %>%
  mutate(K = paste0("K=", 6))

ggplot(df_k2, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k2$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

ggplot(df_k3, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k3$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")

p1 <- ggplot(df_k4, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k4$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")
p1

p2 <- ggplot(df_k5, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k5$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")
p2

p3 <- ggplot(df_k6, aes(x = id, y=Q, fill=reorder(Source, -Q))) +
  geom_col(position = position_stack(), width=1) +
  scale_fill_manual(values = c("lightyellow","#a6a0ff","blue","#008080","yellow", "orange"))+
  #                   labels = c("Native", "AQ_2","AQ_1A","AQ_4", "AQ_3")) +
  geom_vline(xintercept = cumulative+0.5, size=0.25, linetype="dashed") +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("", ""),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="", title=unique(df_k6$K))+
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        plot.title = element_text(size = 16))+
  guides(fill = "none")
p3

# get vline intercepts
table(str_n1601$Grp)
counts <- as.vector(table(str_n1601$Grp))
# cumulative sum
cumulative <- cumsum(counts)
cumulative

