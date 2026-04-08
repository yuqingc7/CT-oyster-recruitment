setwd("/local/workdir/yc2644/CV_CT_array/structure_K/test_recomsim/old_structure_way")

rmsd_all <- bind_rows(
  rmsd_results_a %>% mutate(dataset = "4a"),
  rmsd_results_b %>% mutate(dataset = "4b"),
  rmsd_results_c %>% mutate(dataset = "S5b"),
  rmsd_results_d %>% mutate(dataset = "S5a"),
  rmsd_results_e %>% mutate(dataset = "4c"),
  rmsd_results_f %>% mutate(dataset = "4d")
)
  
rmsd_all$pop <- factor(rmsd_all$pop, levels = unique(rmsd_results_f$pop))

ggplot(rmsd_all, aes(x = pop, y = rmsd, group = dataset, color = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Population", y = "RMSD", color = "Figure") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )
ggplot(rmsd_all, aes(x = pop, y = rmsd, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#457b9d", "#e63946", "#2a9d8f", 
                               "#f4a261", "#8d99ae", "#264653")) +
  labs(x = "Population", y = "RMSD", fill = "Method") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# top 5 905 dapc2 fst2 no centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("test_selected5_905_nocentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V6 = as.numeric(V6),Q = 1-V6) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      str_detect(ind, "WildUnadmix") ~ "WildUnadmixed",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed"))) %>% 
  select(-V1)

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

rmsd_results_d <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_d

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

# top 5 905 dapc2 fst2 only centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("test_selected5_905_onlycentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V4 = as.numeric(V4),Q = V4) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      str_detect(ind, "WildUnadmix") ~ "WildUnadmixed",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  filter(pop %in% c("HR","AQ_1","AQ_2","AQ_3",
                    "AQ_4","AQ_5","AQ_6",
                    "F1HYB","B2Wild","B3Wild",
                    "B4Wild","B5Wild", "WildUnadmixed")) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed"))) %>% 
  select(-V1)

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

rmsd_results_e <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_e

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
  ) +
  ylim(0,1)

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

# 10k random no centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("10kSNP.test_fp_noinvers.fixed.sorted_nocentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V6 = as.numeric(V6),Q = 1-V6) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      str_detect(ind, "WildUnadmix") ~ "WildUnadmixed",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed"))) %>% 
  select(-V1)

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

rmsd_results_a <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_a

ggplot(all_query_Q, aes(x = pop, y = Q)) +
  geom_jitter(shape=1,width = 0.25, size = 2.5, alpha = 0.5, color = "#1d3557") +
  geom_boxplot(color = "black", fill=NA,outliers = FALSE) +
  geom_point(data = expected_vals, aes(x = pop, y = true_Q),
             color = "red", size = 3, alpha=0.5) +
  labs(x = "Population", y = "Introgression Level Q") +
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
# 335 422

# 10k random only centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("10kSNP.test_fp_noinvers.fixed.sorted_onlycentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V4 = as.numeric(V4),Q = V4) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      str_detect(ind, "WildUnadmix") ~ "WildUnadmixed",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  filter(pop %in% c("HR","AQ_1","AQ_2","AQ_3",
                    "AQ_4","AQ_5","AQ_6",
                    "F1HYB","B2Wild","B3Wild",
                    "B4Wild","B5Wild", "WildUnadmixed")) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmixed"))) %>% 
  select(-V1)

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

rmsd_results_b <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results_b

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
  ) +
  ylim(0,1)

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

# 905 random no centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("905SNP.test_fp_noinvers.fixed.sorted_nocentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V4 = as.numeric(V4),Q = 1-V4) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmix"))) %>% 
  select(-V1)

table(all_query_Q$pop)

expected_vals <- tibble(
  pop = factor(c("F1HYB", "B2Wild", "B3Wild", "B4Wild", "B5Wild", "WildUnadmix"),
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

rmsd_results <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results

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

# 905 random only centroids ------------------------------------------------------------------
# Read raw lines
lines <- read_lines("905SNP.test_fp_noinvers.fixed.sorted_onlycentroids_q")
# Split each line by whitespace or tab, remove empty strings
split_lines <- strsplit(lines, "\\s+")  # splits on one or more spaces
split_lines <- lapply(split_lines, function(x) x[x != ""])  # remove empty entries
all_query_Q <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
  select(-V2) %>% 
  mutate(ind = V1, V4 = as.numeric(V4),Q = V4) %>%
  mutate(
    pop = case_when(
      str_detect(ind, "pop1") ~ "HR",
      str_detect(ind, "Cv5786") ~ "AQ_2",
      str_detect(ind, "FI1012|FIS") ~ "AQ_1",
      str_detect(ind, "NEH1|NEH2") ~ "AQ_3",
      str_detect(ind, "UMFS") ~ "AQ_4",
      str_detect(ind, "KRB16") ~ "AQ_5",
      str_detect(ind, "KRP13") ~ "AQ_6",
      str_detect(ind, "F1HYB") ~ "F1HYB",
      str_detect(ind, "B2Wild") ~ "B2Wild",
      str_detect(ind, "B3Wild") ~ "B3Wild",
      str_detect(ind, "B4Wild") ~ "B4Wild",
      str_detect(ind, "B5Wild") ~ "B5Wild",
      TRUE ~ ind
    ),
    grp = case_when(
      str_detect(pop, "^AQ") ~ "AQ",
      TRUE ~ pop
    )
  ) %>% 
  filter(pop %in% c("HR","AQ_1","AQ_2","AQ_3",
                    "AQ_4","AQ_5","AQ_6",
                    "F1HYB","B2Wild","B3Wild",
                    "B4Wild","B5Wild", "WildUnadmix")) %>% 
  mutate(pop=factor(pop, levels=c("HR","AQ_1","AQ_2","AQ_3",
                                  "AQ_4","AQ_5","AQ_6",
                                  "F1HYB","B2Wild","B3Wild",
                                  "B4Wild","B5Wild", "WildUnadmix"))) %>% 
  select(-V1)

table(all_query_Q$pop)

expected_vals <- tibble(
  pop = factor(c("F1HYB", "B2Wild", "B3Wild", "B4Wild", "B5Wild", "WildUnadmix"),
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

rmsd_results <- all_query_Q %>%
  inner_join(expected_vals, by = "pop") %>%
  mutate(sq_diff = (Q - true_Q)^2) %>%
  summarise(rmsd = sqrt(mean(sq_diff)), .by = pop) %>%
  arrange(pop)

rmsd_results

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
  ) +
  ylim(0,1)

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
