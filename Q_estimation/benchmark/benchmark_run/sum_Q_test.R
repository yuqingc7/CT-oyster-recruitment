setwd("/local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/")

# top 5 905 dapc2 fst2 ------------------------------------------------------------------
setwd("/local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/")

all_query_Q <- read_tsv("./top5_dapc_fst/output/rep_0/all_query_Q_values_with_pop.tsv",
                        col_names = FALSE) %>%
  rename(pop = X1, ind = X2, Q = X3) %>%
  mutate(
    pop = case_when(
      str_detect(pop, "HH1013a|HH1013s|HRP22") ~ "HR",
      str_detect(pop, "Cv5786") ~ "AQ_2",
      str_detect(pop, "FI1012|FIS1013a") ~ "AQ_1",
      str_detect(pop, "NEH1|NEH2") ~ "AQ_3",
      str_detect(pop, "UMFS") ~ "AQ_4",
      str_detect(pop, "KRB16") ~ "AQ_5",
      str_detect(pop, "KRP13") ~ "AQ_6",
      TRUE ~ pop
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed")))

table(all_query_Q$pop)

expected_vals <- tibble(
  pop = factor(c("F1HYB", "B2Wild", "B3Wild", "B4Wild", "B5Wild", "WildUnadmixed"),
               levels = levels(all_query_Q$pop)),
  true_Q = c(0.5, 0.25, 0.125, 0.0625, 0.03125, 0)
)

# Compute mean squared difference (MSD)
msd_results <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results

rmsd_results_f <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_f

ggplot(all_query_Q, aes(x = pop, y = Q)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

ggplot(subset(all_query_Q, !(grp %in% c("AQ", "HR"))), aes(x = pop, y = Q)) +
  geom_jitter(shape=1,width = 0.25, size = 2.5, alpha = 0.5, color = "#1d3557") +
  geom_boxplot(color = "black", fill=NA,outliers = FALSE) +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.75, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  ) +
  ylim(0,1)

all_query_Q %>% filter(grp %in% c("HR", "AQ")) %>% 
  group_by(grp) %>% 
  summarize(mean=mean(Q))
# # A tibble: 2 × 2
# grp      mean
# <chr>   <dbl>
#   1 AQ    0.937  
# 2 HR    0.00115

PD <- 0.937
PW <- 0.00115
all_query_Q_rel <- all_query_Q %>% filter(!grp %in% c("HR", "AQ")) %>% 
  mutate(Q_relative=(Q-PW)/(PD-PW)) 
ggplot(all_query_Q_rel, aes(x = pop, y = Q_relative)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Relative Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )


# Compute mean squared difference (MSD)
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff)) 
# All MSD = 0.00107
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)))
# All RMSD = 0.0328

msd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results_rel

# Compute RMSD by population
rmsd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_rel

# Proportion of admixed inds in the pop
# true
table(all_query_Q_rel$pop)
100/120
# 0.8333333
dim(subset(all_query_Q_rel,Q_relative>0))
98/120
# 0.8166667

# Admixture levels for admixed inds
# true
table(all_query_Q_rel$pop)
(0.5*20+0.25*20+0.125*20+0.0625*20+0.0312*20)/100
# 0.19374
mean(subset(all_query_Q_rel,Q_relative>0)$Q_relative)
# 0.2020165

# random 905 ------------------------------------------------------------------
setwd("/local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/")

all_query_Q <- read_tsv("./random_markers/output/rep_0/all_query_Q_values_with_pop.tsv",
                        col_names = FALSE) %>%
  rename(pop = X1, ind = X2, Q = X3) %>%
  mutate(
    pop = case_when(
      str_detect(pop, "HH1013a|HH1013s|HRP22") ~ "HR",
      str_detect(pop, "Cv5786") ~ "AQ_2",
      str_detect(pop, "FI1012|FIS1013a") ~ "AQ_1",
      str_detect(pop, "NEH1|NEH2") ~ "AQ_3",
      str_detect(pop, "UMFS") ~ "AQ_4",
      str_detect(pop, "KRB16") ~ "AQ_5",
      str_detect(pop, "KRP13") ~ "AQ_6",
      TRUE ~ pop
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed")))

table(all_query_Q$pop)

expected_vals <- tibble(
  pop = factor(c("F1HYB", "B2Wild", "B3Wild", "B4Wild", "B5Wild", "WildUnadmixed"),
               levels = levels(all_query_Q$pop)),
  true_Q = c(0.5, 0.25, 0.125, 0.0625, 0.03125, 0)
)

# Compute mean squared difference (MSD)
msd_results <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results

rmsd_results_c <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_c

ggplot(all_query_Q, aes(x = pop, y = Q)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

ggplot(subset(all_query_Q, !(grp %in% c("AQ", "HR"))), aes(x = pop, y = Q)) +
  geom_jitter(shape=1,width = 0.25, size = 2.5, alpha = 0.5, color = "#1d3557") +
  geom_boxplot(color = "black", fill=NA,outliers = FALSE) +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "", y = "Introgression Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.75, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  ) +
  ylim(0,1)

all_query_Q %>% filter(grp %in% c("HR", "AQ")) %>% 
  group_by(grp) %>% 
  summarize(mean=mean(Q))
# # A tibble: 2 × 2
# grp      mean
# <chr>   <dbl>
#   1 AQ    0.937  
# 2 HR    0.00115

PD <- 0.937
PW <- 0.00115
all_query_Q_rel <- all_query_Q %>% filter(!grp %in% c("HR", "AQ")) %>% 
  mutate(Q_relative=(Q-PW)/(PD-PW)) 
ggplot(all_query_Q_rel, aes(x = pop, y = Q_relative)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Relative Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )


# Compute mean squared difference (MSD)
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff)) 
# All MSD = 0.00107
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)))
# All RMSD = 0.0328

msd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results_rel

# Compute RMSD by population
rmsd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_rel

# Proportion of admixed inds in the pop
# true
table(all_query_Q_rel$pop)
100/120
# 0.8333333
dim(subset(all_query_Q_rel,Q_relative>0))
98/120
# 0.8166667

# Admixture levels for admixed inds
# true
table(all_query_Q_rel$pop)
(0.5*20+0.25*20+0.125*20+0.0625*20+0.0312*20)/100
# 0.19374
mean(subset(all_query_Q_rel,Q_relative>0)$Q_relative)
# 0.2020165


# reps top 5 905 dapc2 fst2------------------------------------------------------------------

all_query_Q <- read_tsv("./top5_dapc_fst/output/all_query_Q_values_with_pop.tsv",
                        col_names = FALSE) %>%
  rename(pop = X1, ind = X2, Rep = X3, Q = X4) %>%
  mutate(
    pop = case_when(
      str_detect(pop, "HH1013a|HH1013s|HRP22") ~ "HR",
      str_detect(pop, "Cv5786") ~ "AQ_2",
      str_detect(pop, "FI1012|FIS1013a") ~ "AQ_1",
      str_detect(pop, "NEH1|NEH2") ~ "AQ_3",
      str_detect(pop, "UMFS") ~ "AQ_4",
      str_detect(pop, "KRB16") ~ "AQ_5",
      str_detect(pop, "KRP13") ~ "AQ_6",
      TRUE ~ pop
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed")))

table(all_query_Q$pop)

expected_vals <- tibble(
  pop = factor(c("F1HYB", "B2Wild", "B3Wild", "B4Wild", "B5Wild", "WildUnadmixed"),
               levels = levels(all_query_Q$pop)),
  true_Q = c(0.5, 0.25, 0.125, 0.0625, 0.03125, 0)
)

# Compute mean squared difference (MSD)
msd_results <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results

ggplot(all_query_Q, aes(x = pop, y = Q)) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  facet_wrap(~Rep, nrow=2) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) + #, color = "#2C3E50"
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )

all_query_Q %>% filter(grp %in% c("HR", "AQ")) %>% 
  group_by(grp, Rep) %>% 
  summarize(mean=mean(Q))

all_query_Q %>% filter(grp %in% c("HR", "AQ")) %>% 
  group_by(grp) %>% 
  summarize(mean=mean(Q))
# grp      mean
# <chr>   <dbl>
# 1 AQ    0.949  
# 2 HR    0.00117

PD <- 0.949
PW <- 0.00117
all_query_Q_rel <- all_query_Q %>% filter(!grp %in% c("HR", "AQ")) %>% 
  mutate(Q_relative=(Q-PW)/(PD-PW)) 
ggplot(all_query_Q_rel, aes(x = pop, y = Q_relative)) +
  facet_wrap(~Rep, nrow=2) +
  geom_boxplot(fill = NA, color = "black", outliers = FALSE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "#2C3E50") +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Relative Admixture Level Q") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold")
  )


# Compute mean squared difference (MSD)
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff)) 
# All MSD = 0.00128
all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)))
# All RMSD = 0.0358

msd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(msd = mean(sq_diff), .by = pop) %>%
  arrange(pop)

msd_results_rel

# Compute RMSD by population
rmsd_results_rel <- all_query_Q_rel %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q_relative - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_rel

# Proportion of admixed inds in the pop
# true
table(all_query_Q_rel$pop)
100/120
# 0.8333333
dim(subset(all_query_Q_rel,Q_relative>0))
471/(5*120)
# 0.785

# Admixture levels for admixed inds
# true
table(all_query_Q_rel$pop)
(0.5*20+0.25*20+0.125*20+0.0625*20+0.0312*20)/100
# 0.19374
mean(subset(all_query_Q_rel,Q_relative>0)$Q_relative)
# 0.204698



