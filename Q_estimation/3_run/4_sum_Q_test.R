setwd("/local/workdir/yc2644/CV_CT_array/structure_K/n1601_by_subset")

library(ggpubr)
library(ggsci)

# top 5 905 dapc2 fst2 ------------------------------------------------------------------
meta <- read_tsv("/local/workdir/yc2644/CV_CT_array/vcf_n1601/pop_grp_info_n1601.tsv") %>% 
  select(-Ind) %>% 
  unique()
all_query_Q <- read_tsv("./output/all_query_Q_values_with_pop.tsv",
                        col_names = FALSE) %>%
  rename(pop = X1, ind = X2, Q = X3) %>%
  mutate(
    Pop = case_when(
      str_detect(pop, "SV") ~ "East",
      str_detect(pop, "HH1013a|HH1013s|HRP22") ~ "Hudson",
      str_detect(pop, "Cv5786") ~ "AQ_2",
      str_detect(pop, "FI1012|FIS1013a") ~ "AQ_1",
      str_detect(pop, "NEH1|NEH2") ~ "AQ_3",
      str_detect(pop, "UMFS") ~ "AQ_4",
      str_detect(pop, "KRB16") ~ "AQ_5",
      str_detect(pop, "KRP13") ~ "AQ_6",
      TRUE ~ pop
    )) %>% 
  left_join(meta) %>% select(-pop) %>% 
  mutate(across(everything(), ~ str_replace_all(., "East", "SV"))) %>% 
  mutate(Q=as.double(Q))

pop_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6",  
               "Hudson","SV",
               "TOC0423","TOC23a","GOO0423","GOO23a","NRH0423","LOP0423","SAB0423","HOR",
               "ASC0423","ASC23a", "LIW1",
               "SPT22s", "STP23a", "MIH0423","MCK0423","QUR23a", "EHF0523","EHF23a",
               "Cv5785","JAC0523","JAC23a","HDC0523","HDC23a",
               "FEC0423","FEC23a","HAM0423","BAK0423","LIW2","YBD")
loc_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6",  
               "Hudson","SV",
               "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
               "SPT","MIH","MCK","QUR", "EHF","Cv5785","JAC","HDC",
               "FEC","HAM","BAK","LIW2","YBD")
grp_order <- c("AQ_1","AQ_2","AQ_3","AQ_4","AQ_5","AQ_6", "Hudson", "SV", "CT")

all_query_Q <- all_query_Q %>%
  mutate(Location=factor(Location, levels = loc_order),
         Grp=factor(Grp, levels = grp_order),
         Pop=factor(Pop, levels = pop_order)) %>%
  arrange(Grp,Pop) 

table(all_query_Q$Pop)
table(all_query_Q$Grp)

all_query_Q <- all_query_Q %>% 
  mutate(
    grp = case_when(
      str_detect(Grp, "^AQ") ~ "AQ",
      TRUE ~ Grp
    )) 

ggplot(all_query_Q, aes(x = Pop, y = Q)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Pind") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

ggplot(all_query_Q, aes(x = Location, y = Q)) +
  geom_boxplot(color = "black", fill="white",outliers = FALSE) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.2, color = "#B31B1B") +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Pind") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )


ggplot(subset(all_query_Q, grp %in% c("AQ", "Hudson")), aes(x = Location, y = Q)) +
  geom_jitter(shape=1,width = 0.25, size = 2.5, alpha = 0.5, color = "#B31B1B") +
  geom_boxplot(color = "black", fill=NA,outliers = FALSE) +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  labs(x = "Reference Population", y = "Pind") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(size = 16),
    #axis.title = element_text()
  ) 

all_query_Q %>%
  filter(grp %in% c("Hudson", "AQ")) %>% 
  group_by(grp) %>% 
  summarize(mean=mean(Q))
# grp       mean
# <chr>    <dbl>
# 1 AQ     0.949  
# 2 Hudson 0.00117

PD <- 0.949
PW <- 0.00117
all_query_Q_rel <- all_query_Q %>% filter(!grp %in% c("Hudson", "AQ")) %>% 
  mutate(Q_relative=(Q-PW)/(PD-PW)) %>% 
  mutate(
    Geo3 = case_when(
      Location %in% c("SV",
                      "TOC","GOO","NRH","LOP","SAB") ~ "West",
      Location %in% c("HOR","ASC","LIW1",
                      "SPT","MIH","MCK","QUR", "EHF") ~ "Central",
      Location %in% c("Cv5785","JAC","HDC",
                      "FEC","HAM","BAK","LIW2","YBD") ~ "East",
      TRUE ~ Location
    )) %>%
  mutate(Geo3=factor(Geo3, levels=c("West", "Central", "East"))) %>% 
  mutate(
  Geo2 = case_when(
    Location %in% c("SV",
                    "TOC","GOO","NRH","LOP","SAB","HOR","ASC","LIW1",
                    "SPT","MIH","MCK","QUR", "EHF") ~ "West",
    Location %in% c("Cv5785","JAC","HDC",
                    "FEC","HAM","BAK","LIW2","YBD") ~ "East",
    TRUE ~ Location
  )) %>% 
  mutate(Geo2=factor(Geo2, levels=c("West", "East")))
# ggplot(all_query_Q_rel, aes(x = Pop, y = Q_relative)) +
#   geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
#   geom_jitter(width = 0.2, size = 1.5, alpha = 0.3, color = "#1d3557") +
#   # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
#   #            color = "red", size = 3, alpha=0.5) +
#   labs(x = "Population", y = "Introgression Level Q") +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     plot.title = element_text(face = "bold", size = 16),
#     axis.title = element_text(face = "bold")
#   )

wilcox.test(Q_relative ~ Geo2, data = all_query_Q_rel)
# Wilcoxon rank sum test with continuity correction
# 
# data:  Q_relative by Geo2
# W = 164086, p-value = 9.145e-05
# alternative hypothesis: true location shift is not equal to 0

kruskal.test(Q_relative ~ Location, data = all_query_Q_rel)
# Kruskal-Wallis chi-squared = 181.07, df = 21, p-value < 2.2e-16
# at least one location differs

compare_means(
  formula = Q_relative ~ Geo3,
  #group.by = "Location",
  data = all_query_Q_rel,
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

model_glmm <- glm(Q_relative ~ Location,
                   data = all_query_Q_rel, family = Gamma(link = "log"))
summary(model_glmm)


ggplot(all_query_Q_rel, aes(x = Location, y = Q_relative)) +
  geom_jitter(shape=1,width = 0.25, size = 2.5, alpha = 0.5, color = "#1d3557") +
  geom_boxplot(color = "black", fill=NA,outliers = FALSE) +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  ) 

ggplot(all_query_Q_rel, aes(x = Geo3, y = Q_relative)) +
  geom_boxplot(color = "black", fill="white",outliers = FALSE) +
  geom_jitter(shape=1,width = 0.2, size = 2.5, alpha = 0.5, color = "#1d3557") +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  geom_pwc(method = "wilcox.test",
           label = "{p.adj.format}{p.adj.signif}",
           p.adjust.method = "bonferroni", p.adjust.by = "panel",
           #hide.ns = TRUE
  ) +
  labs(x = "Population", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

compare_means(
  formula = Q_relative ~ Geo2,
  #group.by = "Location",
  data = all_query_Q_rel,
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

all_query_Q_rel %>% 
  group_by(Geo3) %>% 
  summarize(median=median(Q_relative),
            mean=mean(Q_relative))


library(smplot2)
ggplot(all_query_Q_rel, aes(x = Geo3, y = Q_relative, fill=Geo3)) +
  #geom_point(shape=1, size = 2.5, alpha = 0.5, color = "#1d3557") +
  #geom_violin(color = "black", aes(fill=Geo2)) +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),
               point.params = list(
                 size = 2, shape = 21, color = "transparent")) +
  geom_pwc(method = "wilcox.test",
           label = "{p.adj.format}{p.adj.signif}",
           p.adjust.method = "bonferroni", p.adjust.by = "panel",
           hide.ns = TRUE
  ) +
  scale_fill_manual(values = c("#da73e6","lightgreen","#ffe764")) +
  #scale_fill_brewer(palette = "Pastel1") +
  labs(x = "", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
  )

ggplot(all_query_Q_rel, aes(x = Geo2, y = Q_relative, fill=Geo2)) +
  #geom_point(shape=1, size = 2.5, alpha = 0.5, color = "#1d3557") +
  #geom_violin(color = "black", aes(fill=Geo2)) +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),
               point.params = list(
                 size = 2, shape = 21, color = "transparent")) +
  geom_pwc(method = "wilcox.test",
           label = "{p.adj.format}{p.adj.signif}",
           p.adjust.method = "bonferroni", p.adjust.by = "panel",
           #hide.ns = TRUE
  ) +
  scale_fill_manual(values = c("#ffe764","#911eb4")) +
  #scale_fill_brewer(palette = "Pastel1") +
  labs(x = "", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
  )

all_query_Q_rel %>% 
  group_by(Geo2) %>% 
  summarise(median=median(Q_relative))

all_query_Q_rel %>% 
  filter(Q_relative<0.9) %>% 
  group_by(Geo2) %>% 
  summarise(median=median(Q_relative))


all_query_Q_rel %>% 
  filter(Q_relative<0.9) %>% 
  ggplot(aes(x = Geo2, y = Q_relative, fill=Geo2)) +
  #geom_point(shape=1, size = 2.5, alpha = 0.5, color = "#1d3557") +
  #geom_violin(color = "black", aes(fill=Geo2)) +
  # geom_point(data = expected_vals, aes(x = pop, y = true_Q),
  #            color = "red", size = 3, alpha=0.5) +
  sm_raincloud(boxplot.params = list(outlier.shape = NA),
               point.params = list(
                 size = 2, shape = 21, color = "transparent")) +
  geom_pwc(method = "wilcox.test",
           label = "{p.adj.format}{p.adj.signif}",
           p.adjust.method = "bonferroni", p.adjust.by = "panel",
           #hide.ns = TRUE
  ) +
  scale_fill_manual(values = c("#ffe764","#911eb4")) +
  #scale_fill_brewer(palette = "Pastel1") +
  labs(x = "", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
  )
# 
# "#7f404a" "#5b4080" "#408073" "#8c994d" "#cc9666" "#cc1489" "#1262b3" "#cc3d3d" "#da73e6" 
# "#66b1cc" "#0f993d"
# "#7f5d0d" "#7b3dcc" "#45e0e6" "#63e617" "#e57717" "#c9b9c6" "#ffe764" "#ffb359" "#9ee1a8"

#write_tsv(all_query_Q_rel, "all_query_Q_rel.tsv")

all_query_Q_rel %>% 
  ggplot(aes(x = Q_relative,fill=Geo2)) +
  facet_wrap(~Geo2, nrow=2) +
  geom_histogram(aes(y = stat(count / sum(count))),binwidth=0.1) +
  #geom_histogram(binwidth=0.1) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
  ) + coord_cartesian(xlim=c(0,1))

# Proportion of admixed inds in the pop
dim(all_query_Q_rel)
dim(subset(all_query_Q_rel,Q_relative>0))
1317/1317
# 1

# Admixture levels for admixed inds
mean(subset(all_query_Q_rel,Q_relative>0)$Q_relative)
# 0.2416288
range(subset(all_query_Q_rel,Q_relative>0)$Q_relative)
median(subset(all_query_Q_rel,Q_relative>0)$Q_relative)

# adult vs spat ------------------------------------------------------------

loc_order <- c("TOC","GOO","ASC","SPT","EHF","JAC","HDC","FEC")
pop_order <- c("TOC0423","TOC23a","GOO0423","GOO23a",
               "ASC0423","ASC23a", "SPT22s","STP23a","EHF0523","EHF23a",
               "JAC0523","JAC23a", "HDC0523", "HDC23a","FEC0423", "FEC23a")

CT_pair <- all_query_Q_rel %>% 
  filter(Location %in% c("ASC","EHF","FEC","GOO",
                         "HDC","JAC","SPT","TOC")) %>% 
  mutate(Location=factor(Location, levels = loc_order),
         Pop=factor(Pop,levels=pop_order)) %>%
  mutate(Type = recode(Type, "CT_adult" = "Adult", "CT_spat" = "Spat"),
         Type = factor(Type, levels = c("Spat", "Adult")))

library(lme4)
library(lmerTest)  # for p-values

# Each Location has a fixed, separate effect on AQ
# Estimates effect of each location
summary(lm(Q_relative ~ Type + Location, data = CT_pair)) 

# Location effects are random samples from a larger population of locations.
# Estimates variance among locations
model <- lmer(Q_relative ~ Type + (1 | Location), data = CT_pair)
summary(model)

res <- residuals(model)
hist(res, breaks = 30, main = "Histogram of residuals", xlab = "Residuals")
qqnorm(res)
qqline(res, col = "red")
shapiro.test(res)
# residuals are not normal

# beta regression?
# Q is bounded (0–1). lmer assumes Gaussian residuals, 
# which may not hold if Q is near 0 or 1 or highly skewed.
hist(CT_pair$Q_relative)
shapiro.test(CT_pair$Q_relative)
# actually now Q is pretty normally distributed

# Use a generalized model
# If your response is bounded (e.g., proportion, rate, or ratio), 
# a generalized linear mixed model (GLMM) might be more appropriate:
library(lme4)
model_glmm1 <- glm(Q_relative ~ Type + Location,
                    data = CT_pair, family = Gamma(link = "log"))
summary(model_glmm1)

model_glmm2 <- glmer(Q_relative ~ Type + (1 | Location),
                    data = CT_pair, family = Gamma(link = "log"))
summary(model_glmm2)

library(lattice)  # for better plotting

# Extract residuals
res <- residuals(model_glmm, type = "pearson")  # or "deviance"
fitted_vals <- fitted(model_glmm)

hist(res, breaks = 30, main = "Histogram of Pearson Residuals",
     xlab = "Residuals")

qqnorm(res)
qqline(res, col = "red")

plot(fitted_vals, res,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")


library(ggpubr)
library(ggsci)
# ggplot(CT_pair, aes(x = Type, y = Q_relative, fill = Type)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
#   geom_jitter(position = position_jitter(width = 0.15), 
#               size = 2, alpha = 0.5, shape = 21, color = "black") +
#   stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = T,
#                      size = 8,vjust = -0.2,label.y.npc = "top",fontface = "bold") +
#   scale_fill_brewer(palette = "Set1") +
#   #coord_cartesian(ylim = c(0, 0.5)) +
#   facet_wrap(~ Location, nrow=1) +
#   theme_minimal(base_size = 16) +
#   labs(title = "",
#        x = "", y = "Q", fill = "Type") +
#   theme(
#     strip.text = element_text(face = "bold", size = 14),
#     plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
#     legend.position = "none"
#   )
# 
# ggplot(CT_pair, aes(x = Type, y = Q_relative, fill = Type)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
#   geom_jitter(position = position_jitter(width = 0.15), 
#               size = 2, alpha = 0.5, shape = 21, color = "black") +
#   stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = T,
#                      size = 8,vjust = -0.2,label.y.npc = "top",fontface = "bold") +
#   scale_fill_brewer(palette = "Set1") +
#   #coord_cartesian(ylim = c(0, 0.5)) +
#   #facet_wrap(~ Location, nrow=1) +
#   theme_minimal(base_size = 16) +
#   labs(title = "",
#        x = "", y = "Q", fill = "Type") +
#   theme(
#     strip.text = element_text(face = "bold", size = 14),
#     plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
#     legend.position = "none"
#   )

compare_means(
  formula = Q_relative ~ Type,
  group.by = "Location",
  data = CT_pair,
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
)

a <- ggplot(CT_pair, 
       aes(x = Type, y = Q_relative)) +
  geom_jitter(aes(shape = Type), 
              position = position_jitter(width = 0.3), 
              size = 2.5, alpha = 0.6, fill = "blue", color="white") +
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.6, fill = "white") +
  # stat_compare_means(method = "wilcox.test",
  #                    label = "p.signif",
  #                    hide.ns = TRUE,
  #                    size = 8,
  #                    vjust = -0.2,
  #                    label.y.npc = 0.9,
  #                    fontface = "bold") +
  geom_pwc(method = "wilcox.test",
           label = "{p.adj.format}{p.adj.signif}",
           p.adjust.method = "bonferroni", p.adjust.by = "panel",
           hide.ns = TRUE
           ) +
  scale_shape_manual(values = c("Spat" = 24, "Adult" = 21)) + # 24=triangle, 21=circle
  coord_cartesian(ylim = c(0, 0.65)) +
  facet_wrap(~ Location, nrow = 1) +
  labs(x = "", y = "Introgression Level Q", shape = "Type") +
  # theme(
  #   strip.text = element_text(face = "bold", size = 12),
  #   plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
  #   axis.text.x = element_text(size = 10),
  #   legend.position = "none"
  # ) 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggadjust_pvalue(a,
  p.adjust.method = "bonferroni",
  label = "{p.adj.format}{p.adj.signif}", hide.ns = TRUE
) 

compare_means(
  formula = Q_relative ~ Type,
  #group.by = "Type",
  data = CT_pair,
  method = "wilcox.test"
)

ggplot(CT_pair, 
       aes(x = Type, y = Q_relative, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(position = position_jitter(width = 0.15), 
              size = 2, alpha = 0.5, shape = 21, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = T,
                     size = 8,vjust = -0.2,label.y.npc = "top",fontface = "bold") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0, 0.5)) +
  #facet_wrap(~ Location, nrow=1) +
  theme_minimal(base_size = 16) +
  labs(title = "",
       x = "", y = "Q", fill = "Type") +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    legend.position = "none"
  )


# compare runs ------------------------------------------------------------


# compare with STRUCTURE naive runs
n1601_k5 <- readQ("/local/workdir/yc2644/CV_CT_array/structure/n1601_cor/str_K5_rep1_f", filetype="auto")$str_K5_rep1_f
grps <- read_tsv("/local/workdir/yc2644/CV_CT_array/vcf_n1601/pop_grp_info_n1601.tsv")
dim(grps)
# combine sample IDs to qlist
str_n1601 <- cbind(grps, n1601_k5) %>% 
  mutate(AQ=1-Cluster4) %>% select(Pop, Ind, Batch, Location, Grp, Type, AQ) %>% 
  mutate(ind = str_replace_all(Ind, "_", "-")) %>% 
  filter(Grp %in% c("CT", "East"))

compare <- left_join(all_query_Q_rel, str_n1601) %>% 
  select(Ind, Pop, Batch, Location, Grp, Type, grp, Q_relative, AQ)
compare %>% 
  ggplot() +
  #facet_wrap(~Pop) +
  geom_point(aes(x=AQ, y=Q_relative), alpha=0.8, shape=21, size=3) + 
  geom_abline(slope=1,intercept=0, color="red")+
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  labs(x="STRUCTURE run of n=1601 individuals (K=5) Pind",
       y="Karlsson Introgression Levels Q", title="1317 Unknown Wild Inds")

ggscatter(compare, x = "AQ", y = "Q_relative",
          color="black", shape = 1, # for points
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, xlim=c(0,1),ylim=c(0,1),
          cor.coef = TRUE, cor.method = "pearson") +
  #geom_vline(xintercept=0.007, color="red", linetype="dashed") +# threshold for admixed with STRUCTURE
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="red") +
  labs(x="STRUCTURE run of n=1601 individuals (K=5) Pind",
       y="Karlsson Introgression Levels Q", title="1317 Unknown Wild Inds")


# ELAI results
tab_graph<- read.table("/local/workdir/yc2644/CV_CT_array/ELAI/06_tracts/Intro_byInd_threshold_0.05.txt", 
                       header=T) %>% 
  select(-Ind, -Pop)
tab_graph
colnames(cbind(str_n1601,tab_graph))
compare_elai <- cbind(str_n1601,tab_graph) %>% 
  right_join(all_query_Q_rel) %>% 
  mutate(m.total_percent_intro=m.total_percent_intro/100) %>% 
  select(Ind, Pop, Batch, Location, Grp, Type, grp, Q, Q_relative, m.total_percent_intro)
compare_elai %>% 
  ggplot() +
  #facet_wrap(~Pop) +
  geom_point(aes(x=m.total_percent_intro, y=Q_relative), alpha=0.8, shape=21, size=3) + 
  geom_abline(slope=1,intercept=0, color="red")+
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1)) +
  labs(x="ELAI Introgression Levels Q",
       y="Karlsson Introgression Levels Q", title="1317 Unknown Wild Inds")

cor.test(compare_elai$Q_relative, compare_elai$m.total_percent_intro, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  compare_elai$Q_relative and compare_elai$m.total_percent_intro
# t = 34.128, df = 1315, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6556008 0.7129694
# sample estimates:
#      cor 
# 0.685347 

hist(compare_elai$Q_relative,col="lightblue")
qqnorm(compare_elai$Q_relative)
qqline(compare_elai$Q_relative, col="red")

ggscatter(compare_elai, x = "m.total_percent_intro", y = "Q_relative",
          color="black", shape = 1, # for points
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, xlim=c(0,1),ylim=c(0,1),
          cor.coef = TRUE, cor.method = "pearson") +
  #geom_vline(xintercept=0.007, color="red", linetype="dashed") +# threshold for admixed with STRUCTURE
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="red") +
  labs(x="ELAI Pind",
       y="Karlsson Introgression Levels Q", title="1317 Unknown Wild Inds")

ggscatter(compare_elai, x = "m.total_percent_intro", y = "Q",
          color="black", shape = 1, # for points
          add = "reg.line",
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE, xlim=c(0,1),ylim=c(0,1),
          cor.coef = TRUE, cor.method = "pearson") +
  #geom_vline(xintercept=0.007, color="red", linetype="dashed") +# threshold for admixed with STRUCTURE
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color="red") +
  labs(x="ELAI Pind",
       y="Karlsson Pind", title="1317 Unknown Wild Inds")

