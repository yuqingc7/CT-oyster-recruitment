setwd("/local/workdir/yc2644/CV_CT_array")

## install.packages("devtools")
#devtools::install_github("GuangchuangYu/scatterpie")
library(scatterpie)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(wesanderson)


# get percent admixture for each site
str_n1367 <- read_tsv("./structure_K/n1601_by_subset/all_query_Q_rel.tsv")

table(str_n1367$Grp)

str_n1367 <- str_n1367 %>%
  mutate(Class = case_when(
    Q_relative > 0  & Q_relative <= 0.25 ~ "Class1",
    Q_relative <= 0.5 & Q_relative > 0.25 ~ "Class2",
    Q_relative > 0.5 ~ "Class3",
    #Q_relative > 0.5 ~ "Class4",
    TRUE ~ "Other"
  ))

unique(str_n1367$Class)

table(str_n1367$Location)

#CT_ER_coord <- read_tsv("CT_ER_coordinates_new.tsv", col_names = F)
CT_ER_coord <- read_tsv("CT_ER_coordinates_jitter_new.tsv", col_names = F) %>% 
  mutate(across(everything(), ~ str_replace_all(., "East", "SV"))) %>% 
  mutate(X2=as.double(X2),X3=as.double(X3))

df <- left_join(str_n1367,CT_ER_coord, join_by(Location == X1)) %>% 
  rename(lat=X2, long=X3)
df[rowSums(is.na(df)) > 0, ]

df_results <- df %>%
  group_by(Location, Class) %>%
  dplyr::summarize(tempCount = n()) 
df_results

df_results_wider <- df_results %>%
  group_by(Location) %>%
  # dplyr::summarize(Count = sum(tempCount),
  #                  Unadmixed = sum(tempCount[Class == "Class1"]),
  #                  `<=12.5 % Admixed` = sum(tempCount[Class == "Class2"]),
  #                  `12.5-25 % Admixed` = sum(tempCount[Class == "Class3"]),
  #                  `25-50 % Admixed` = sum(tempCount[Class == "Class4"]),
  #                  `50-75 % Admixed` = sum(tempCount[Class == "Class5"]),
  #                  `>75 % Admixed` = sum(tempCount[Class == "Class6"]))  %>%
  dplyr::summarize(Count = sum(tempCount),
                   #Unadmixed = sum(tempCount[Class == "Class1"]),
                   `Introgression <=25 %` = sum(tempCount[Class == "Class1"]),
                   `Introgression 25-50 %` = sum(tempCount[Class == "Class2"]),
                   `Introgression >50 %` = sum(tempCount[Class == "Class3"]))  %>%
  ungroup()
# df_results_wider[, 3:8] <- df_results_wider[, 3:8]/df_results_wider$Count
df_results_wider[, 3:5] <- df_results_wider[, 3:5]/df_results_wider$Count
df_results_wider

df_results_wider <- left_join(df_results_wider,CT_ER_coord,join_by(Location == X1)) %>% 
  rename(lat=X2, long=X3)
df_results_wider

colnames(df_results_wider)

# Get high-resolution state boundaries
states_hires <- ne_states(country = "United States of America", 
                          returnclass = "sf")

# Convert to data frame for ggplot
states_df <- fortify(states_hires)

# p <- ggplot(usa, aes(x = long, y = lat, group = group)) +
#   coord_cartesian(xlim = c(-74.05, -72), ylim = c(40.55, 41.5))+
#   geom_polygon(color = "black",fill = "#F0EAD6") +
#   labs(x = "Longitude",y = "Latitude") +
#   theme_minimal()+
#   # Modifying text size for axis titles and ticks
#   theme(axis.title = element_text(size = 16),  # Adjust the font size for axis titles
#         axis.text = element_text(size = 14)    # Adjust the font size for axis ticks
#   )
# #p
# p + geom_scatterpie(data=df_results_wider,
#                     aes(x=long, y=lat,r=0.01*log(Count)),
#                     cols=c("Unadmixed","<=12.5 % Admixed","12.5-25 % Admixed",
#                            "25-50 % Admixed","50-75 % Admixed",">75 % Admixed"),
#                     color=NA,alpha=0.75)+
#   scale_fill_manual(values=col_pal)+
#   guides(fill = guide_legend(title = "Percent Aquaculture-source Admixture Per Sample"))+
#   theme(legend.position = "top")+
#   geom_scatterpie_legend(0.01*log(df_results_wider$Count), n=3,
#                          x=-73.875, y=41.35,labeller=function(x) round(exp(100*x),0))

col_pal <- c("#F1BB7B", "#2E8B57","#FD6467")
range(df_results_wider$Count)
df_results_wider$r = 0.005 * sqrt(df_results_wider$Count)
ggplot() +
  geom_sf(data=states_hires, fill="#a3ada0", color="black", linewidth = 0.4, alpha=0.5) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-74, -72), ylim = c(40.55, 41.45)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  guides(color = "none",fill="none",shape="none") + 
  geom_scatterpie(data=df_results_wider,
                    aes(x=long, y=lat,r=r),
                    # cols=c("Unadmixed","<=12.5 % Admixed","12.5-25 % Admixed",
                    #        "25-50 % Admixed","50-75 % Admixed",">75 % Admixed"),
                    cols=c("Introgression <=25 %",
                           "Introgression 25-50 %","Introgression >50 %"),
                    color=NA,alpha=0.95,sorted_by_radius=T)+
  scale_fill_manual(values = col_pal)+
  #guides(fill = guide_legend(title = "Percent Aquaculture-source Admixture Per Sample"))+
  theme(legend.position = "top")+
  geom_scatterpie_legend(breaks=c(0.005*sqrt(25), 0.005*sqrt(75), 0.005*sqrt(160)),
                         x=-73.75, y=41.35,labeller=function(x) round((x/0.005)^2))

