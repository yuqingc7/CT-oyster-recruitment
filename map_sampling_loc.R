setwd("/local/workdir/yc2644/CV_CT_array")

library(tidyverse)
require(maps)
require(mapdata)

coord <- read_tsv("all_coordinates.tsv", col_names = F)
colnames(coord) <- c("site", "lat", "lon","type")
coord <- coord %>% 
  mutate(type=factor(type, levels=c("spat", "adult", "both")))

w2hr <- map_data("worldHires", "USA")

# ggplot() +
#   geom_polygon(data=w2hr,aes(long,lat,group=group),fill="#a3ada0", color="black",
#                size=0.25, alpha=0.5) +
#   geom_point(data=coord, aes(x = lon, y = lat), size=5,fill="red", color="black",
#              shape=21,alpha=0.5)+
#   labs(x = "Longitude", y = "Latitude") +
#   guides(color="none")+
#   #coord_cartesian(xlim = c(-74.05, -72), ylim = c(40.55, 41.5)) +
#   coord_map(xlim = c(-74, -72), ylim = c(40.55, 41.45), projection = "mercator")+
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14)
#   )

#install.packages(c("rnaturalearth", "rnaturalearthdata", "sf"))
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Get high-resolution state boundaries
states_hires <- ne_states(country = "United States of America", 
                          returnclass = "sf")

# Convert to data frame for ggplot
states_df <- fortify(states_hires)

table(coord$type)
ggplot() +
  geom_sf(data=states_hires, fill="#a3ada0", color="black", linewidth = 0.4, alpha=0.5) +
  geom_point(data=coord, aes(x = lon, y = lat, color=type, shape=type), 
             size=4, stroke=0.7,
             #color="black",
            alpha=0.5) +
  scale_shape_manual(values = c(17,16,8),
                    name="",
                    #  labels = c("CT Adult", "CT Spat", "ER", "HR")
  ) +
  scale_color_manual(values = c("red","red","blue"),
                    name="",
                    #  labels = c("CT Adult", "CT Spat", "ER", "HR")
                    ) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-74, -71.9), ylim = c(40.55, 41.45)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text=element_text(size=14)
  ) 
  #guides(color = "none",fill="none",shape="none")

