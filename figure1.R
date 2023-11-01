# ============================ #
#
# Nephrops RAD Analysis R Script 2023
#
# Figure 1
#
# Data used to create figure:
# ./Data/X.csv
# ./Data/Y.csv
# ./Data/Z.csv
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(tidyverse)
library(sf)
library(terra)
library(ggspatial)
library(RColorBrewer)
library(ggnewscale)
library(ggsflabel)
library(patchwork)
library(jpeg)
library(grid)

# Import site coordinates and convert to sf object
sites <- read_csv("./data/site_coordinates.csv") |>
  st_as_sf(x = _, coords=c("Lon","Lat"), crs=4326) |>
  st_transform(x=_, crs=3035)
med_sites <- filter(sites, Code != "Cly")
clyde <- filter(sites, Code == "Cly")

# Import bathymetry and transform to CRS
bathy <- rast("./data/MS_bathy_5m_lonlat.tif") |>
  terra::project(x = _, "epsg:3035")

# Divide depth into categories
depth_categories <- tibble(bathy=as.numeric(values(bathy)))
depth_categories <- mutate(depth_categories, depth_category = case_when(
  (bathy <= 0 & bathy >= -20) ~ "0-20",
  (bathy <= -20 & bathy >= -40) ~ "20-40",
  (bathy <= -40 & bathy >= -60) ~ "40-60",
  (bathy <= -60 & bathy >= -80) ~ "60-80",
  (bathy <= -80 & bathy >= -100) ~ "80-100",
  (bathy <= -100 & bathy >= -150) ~ "100-150",
  (bathy <= -150 & bathy >= -200) ~ "150-200",
  (bathy <= -200 & bathy >= -1000) ~ "200-1000",
  (bathy < -1000) ~ "> 1000"
))

# Create a new raster with the depth categories
values(bathy) <- depth_categories$depth_category

# Order depth categories
depth_order <- c("0-20","20-40","40-60","60-80","80-100","100-150","150-200","200-1000","> 1000")

# Depth colour palette
blues <- brewer.pal(n_distinct(values(bathy), na.rm = TRUE), "Blues")

#--------------#
# Mediterranean Map
#--------------#

# Mediterranean study area
med_ext_latlon <- ext(7, 26.5, 36, 46)
med_ext_laea <- terra::project(med_ext_latlon, "epsg:4326", "epsg:3035")

# Import world polygons and crop to Mediterranean study area
med_basemap <- st_read("./data/world.gpkg") |>
  st_transform(x=_, crs=3035) |>
  st_crop(x = _, med_ext_laea)

# Import geographical subareas (GSA) polygons and crop to study area
gsa <- st_read("./data/GSAs_simplified/GSAs_simplified.shp") |>
  dplyr::select(.data = _, F_GSA_LIB) |>
  st_transform(x=_, crs=3035)
# plot(gsa)

# Crop bathymetry to study area
med_bathy <- bathy |> crop(x = _, med_ext_laea)
# plot(med_bathy)

# Plot Mediterranean map with bathymetry, GSAs and sampling locations
med_map <- ggplot()+
  layer_spatial(med_bathy, alpha=0.8)+
  scale_fill_manual("Depth (m)", na.translate=FALSE, values=blues, breaks=depth_order)+
  layer_spatial(gsa, fill=NA, aes(linetype=""), linewidth=1, colour="black")+
  scale_linetype_manual("GSA boundaries", values=1)+
  layer_spatial(med_basemap, colour="black", linewidth=0.3)+
  new_scale_fill()+
  layer_spatial(med_sites, aes(fill=GSA), shape=21, colour="black", size=4)+
  geom_sf_label_repel(data=med_sites, aes(label=Code), size = 4,
  force = 30, nudge_x = -2, nudge_y = 0, seed = 10)+
  scale_fill_manual("Sampling location",
                    values=c("yellow","green","#26F7FD","pink"),
                    labels=c("GSA 17","GSA 18","GSA 22","GSA 9"))+
  coord_sf(xlim=c(med_ext_laea[1], med_ext_laea[2]),
           ylim=c(med_ext_laea[3], med_ext_laea[4]),
           expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  annotation_north_arrow(
    data = med_basemap,
    which_north = "true",
    location = "tl",
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm"),
    # pad_y = unit(0.8, "cm"),
    style = north_arrow_orienteering(text_size = 5)
  )+
  # annotation_scale(
  #   data = med_basemap,
  #   location = "tl",
  #   width_hint = 0.2,
  #   bar_cols = c("black","white"),
  #   # height = unit(0.5, "cm"),
  #   text_cex = 1
  # )+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth=1)
  )
# med_map

#--------------#
# Atlantic Map (Firth of Clyde)
#--------------#

# Atlantic study area
atl_ext_latlon <- ext(-9, 28, 32, 59)
atl_ext_laea <- terra::project(atl_ext_latlon, "epsg:4326", "epsg:3035")

# Import world polygons and crop to Atlantic study area
atl_basemap <- st_read("./data/world.gpkg") |>
  st_transform(x=_, crs=3035) |>
  st_crop(x = _, atl_ext_laea)

# Crop bathymetry to study area
atl_bathy <- bathy |> crop(x = _, atl_ext_laea)
# plot(atl_bathy)

# Plot Atlantic map with Firth of Clyde sampling location
atl_map <- ggplot()+
  layer_spatial(atl_bathy, alpha=0.8)+
  scale_fill_manual("Depth (m)", na.translate=FALSE, values=blues, breaks=depth_order)+
  layer_spatial(atl_basemap, colour="black", linewidth=0.3)+
  layer_spatial(clyde, shape=21, fill="white", colour="black", size=4)+
  geom_sf_label_repel(data=clyde, aes(label=Code), size = 4.5,
                      force = 10, nudge_x = 5, seed = 10)+
  geom_rect(aes(xmin=med_ext_laea[1], xmax=med_ext_laea[2],
                ymin=med_ext_laea[3], ymax=med_ext_laea[4]),
            fill=NA, colour="black", linetype=5, linewidth=1)+
  coord_sf(xlim=c(atl_ext_laea[1], atl_ext_laea[2]),
           ylim=c(atl_ext_laea[3], atl_ext_laea[4]),
           expand=FALSE)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth=1),
    legend.position = "none"
  )
  
#--------------#
# Figure 1A: Inset Map
#--------------#

# Superimpose the Atlantic map on the Mediterranean map
fig1A <- med_map+ inset_element(atl_map, 0.6, 0.6, 1, 0.99, align_to = "panel")
# fig1A 

#--------------#
# Figure 1B: Nephrops image
#--------------#

# Import jpeg and convert to ggplot object
nephrops_jpeg = readJPEG("./img/Langoustine_dans_son_terrier.jpg") |>
  rasterGrob(image = _, interpolate = TRUE)
fig1B = ggplot()+
  annotation_custom(nephrops_jpeg, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_void()
# fig1B

#--------------#
# Figure 1C: Violin plot of Nephrops sizes
#--------------#

# Import Nephrops carapace length data
# carapace_df <- read_csv("./data/growth_data_rad_samples_only.csv")
carapace_df <- read_csv("./data/growth_data_all_samples_collected.csv")

# Remove missing data
# carapace_df <- drop_na(carapace_df)

# Reorder sites
site_order <- c("Cly","17I","18II","Anc","Cgg","Pom1","Pom2","Pom3")
carapace_df$Site <- factor(carapace_df$Site, levels = site_order)

# Violin plot grouped by sex (male or female)
fig1C <- ggplot(data=carapace_df, aes(x=Site, y=Carapace_length_mm))+
  # geom_violin(aes(fill=Sex), position = position_dodge(0.7))+
  geom_boxplot(aes(fill=Sex), position = position_dodge(0.7))+
  scale_fill_manual(values= c("#dd1c77","royalblue"),
                    labels= c("Female","Male"))+
  ylab("Carapace length (mm)\n")+
  ggtitle("Carapace length variation")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black")
  )
# fig1C

# ANOVA of carapace lengths in the Adriatic Sea
anova <- carapace_df |>
  filter(.data = _, Sea == "Adriatic") |>
  mutate(.data = _, Pomo = as.factor(Pomo)) |>
  lm(Carapace_length_mm ~ Pomo, data = _)

# ANOVA of carapace lengths in the Adriatic Sea with sex as fixed effect
anova <- carapace_df |>
  filter(.data = _, Sea == "Adriatic") |>
  mutate(.data = _, Pomo = as.factor(Pomo)) |>
  lm(Carapace_length_mm ~ Pomo + Sex, data = _)

# Summary of ANOVA
summary(anova)


#--------------#
# Figure 1 Export
#--------------#

# Layout design
# layout <- "
#   AAAABB
#   AAAABB
#   AAAABB
#   AAAACC
#   AAAACC
#   AAAACC
# "
layout <- "
  AAAA
  AAAA
  BBCC
"

# Plot layout
plt_list = list(
  wrap_elements(fig1A)+ labs(tag = "A"),
  fig1B+ labs(tag = "B"),
  fig1C+ labs(tag = "C")
)
figure1 <- wrap_plots(plt_list, design = layout)

# Export
ggsave(plot = figure1, filename = "figure1.jpg", width = 10, height = 9, dpi = 600)
