# ============================ #
#
# Nephrops RAD Analysis R Script 2023
#
# SNP Population Structure
#
# Data source:
# ./Data/nephrops_snps.RData
# Genind file output from snp_qc_filtering.R script.
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(adegenet)
library(poppr)
library(hierfstat)
library(mmod)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(scales)
library(dartR)
library(SNPRelate)
library(stringr)
library(rlang)
library(tidyr)
library(dplyr)
library(patchwork)

# Source R functions
source("./utils.R")

# Import QC SNPs
load("./data/nephrops_snps.RData"); nephrops_gl

# Import RAD samples metadata
rad_samples_meta <- read.csv("./data/growth_data_rad_samples_only.csv")

# Assess sex ratio per site
rad_samples_meta |>
  filter(.data = _, Ind_ID %in% indNames(nephrops_gl)) |>
  count(x = _, Site, Sex)

#--------------#
# FST ####
#--------------#

# Compute pairwise Fst
Fst <- genet.dist(gl2gi(nephrops_gl), method = "WC84")
# D <- pairwise_D(gl2gi(nephrops_gl))

# Convert minus values to zero
Fst[Fst < 0] <- 0

# Desired order of labels
lab_order <- c("Cly","9I","Aeg","18II","17I","Anc","Cgg","Pom1","Pom2","Pom3")

# Change order of rows and cols
fst.mat <- as.matrix(Fst)
fst.mat1 <- fst.mat[lab_order, ]
fst.mat2 <- fst.mat1[, lab_order]

# Create a data.frame
ind <- which(upper.tri(fst.mat2), arr.ind = TRUE)
fst.df <- data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                     Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                     Fst = fst.mat2[ ind ])

# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 <- factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 <- factor(fst.df$Site2, levels = unique(fst.df$Site2))

# Round digits
fst.df$Fst <- round(fst.df$Fst, digits = 3)

# Fst italic label
fst.label <- expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid <- max(fst.df$Fst) / 2

# Plot heatmap
fst_heat <- ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="white", size = 2)+
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 9, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.2,0.8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 8),
  )
fst_heat

#--------------#
# Perform PCA with Clyde Sea ####
#--------------#

# Impute missing data then run PCA
pca1 <- nephrops_gl |>
  gl.impute(x = _, method = "neighbour") |>
  # tab(x = _, NA.method = "mean") |>
  glPca(x = _, scale = TRUE, center = TRUE, nf = 3)

# Extract percent of genetic variance is explained by each axis
(percent <- pca1$eig/sum(pca1$eig)*100)

# PCA colours
cols = c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#d6604d"
         ,"#2171b5","#2171b5","#2171b5")

# Plot PCA
(pca_plt <- plot_pca(pca1$scores, percent, nephrops_gl, rad_samples_meta, by = "site", cols = cols))
# (pca_plt <- plot_pca(pca1$scores, percent, nephrops_gl, rad_samples_meta, by = "sex", axes = c(1,2)))

#--------------#
# Perform Snapclust with Clyde Sea ####
#--------------#

# Assess number of clusters using AICc
AICc <- snapclust.choose.k(gl2gi(nephrops_gl), max = nPop(nephrops_gl), IC = "AICc")

# Plot AICc
AICc_plt <- ggplot(data = data.frame(AICc = AICc, K = 1:length(AICc)))+
  geom_point(aes(x=K, y=AICc), shape = 19, size = 2.5, colour = "black")+
  geom_line(aes(x=K, y=AICc))+
  scale_x_continuous(breaks = 1:length(AICc))+
  ylab("AICc")+
  xlab("K")+
  # ggtitle("Assess number of clusters using AICc")+
  theme(
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(face = "italic"),
    plot.title = element_text(hjust=0.5, size=15) 
  )
AICc_plt

# K
K <- 2
  
# Impute missing data then run find.clusters()
grp <- nephrops_gl |>
  gl.impute(x = _, method = "neighbour") |>
  find.clusters(x = _, n.clust = K, n.pca = 200)
  
# Run snapclust
snapclust_res <- nephrops_gl |>
  gl.impute(x = _, method = "neighbour") |>
  gl2gi(x = _) |>
  snapclust(x = _, k = K, pop.ini = grp$grp)
# compoplot(snapclust_res)

# Export membership probabilities for mapmixture
mem_prob <- tibble(
  Site = nephrops_gl$pop,
  Ind = names(snapclust_res$group),
  `Cluster 1` = as.data.frame(snapclust_res$proba)$`1`,
  `Cluster 2` = as.data.frame(snapclust_res$proba)$`2`
)
mem_prob
# write.table(mem_prob, file = "./data/snapclust_membership_prob.txt", row.names = FALSE)

#--------------#
# Figure 2 ####
#--------------#

# Read in mapmixture coordinates file
coords <- read.csv("data/mapmixture_coordinates.csv")

# Run mapmixture
library(mapmixture)
mapmixture_plt <- mapmixture(
  admixture_df = mem_prob,
  coords_df = coords,
  cluster_names = c("Cluster 1", "Cluster 2"),
  cluster_cols = c("#2166AC", "#D6604D"),
  boundary = c(xmin=-5.3, xmax=26.5, ymin=34, ymax=58),
  crs = 3035,
  arrow_position = "tl",
  arrow_size = 1.5,
  scalebar_position = "bl",
  scalebar_size = 1.3,
  axis_text_size = 10,
  axis_title_size = 12
)+
  theme(
    legend.text = element_text(size = 12),
    legend.position = "right",
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))

# Layout design
layout <- "
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD#
"

# Plot layout
plt_list = list(
  fst_heat+ labs(tag = "A"),
  pca_plt+ labs(tag = "B")+ theme(axis.title = element_text(size=10)),
  AICc_plt+ labs(tag = "C"),
  mapmixture_plt+ labs(tag = "D")
)
figure2 <- wrap_plots(plt_list, design = layout)

# Export
ggsave(plot = figure2, filename = "figure2.pdf", width = 12, height = 10)

#--------------#
# Perform PCA / DAPC with Med Sites Only ####
#--------------#

# Data set for PCA (remove Clyde Sea and subsequent monomorphic SNPs)
nephrops_med_gl <- nephrops_gl |>
  gl.drop.pop(x = _, pop.list = c("Cly"), mono.rm = TRUE)

# PCA colours
cols = c("#FFC000","#7fcdbb","#bf812d","#dfc27d","#DFFF00","#FDDA0D","yellow","yellow","yellow")

# Impute missing data then run PCA
pca2 <- nephrops_med_gl |>
  gl.impute(x = _, method = "neighbour") |>
  # tab(x = _, NA.method = "mean") |>
  glPca(x = _, scale = TRUE, center = TRUE, nf = 3)

# Extract percent of genetic variance is explained by each axis
(percent <- pca2$eig/sum(pca2$eig)*100)

# Plot PCA
(pca2_plt <- plot_pca(pca2$scores, percent, nephrops_med_gl, rad_samples_meta,
                      by = "site", cols = cols,
                      title = "Principal component analysis: Mediterranean sites"))
# plot_pca(pca2$scores, percent, nephrops_med_gl, rad_samples_meta, by = "sex")

# # Impute missing data then run DAPC
# nephrops_med_gl <- gl.impute(x = nephrops_med_gl, method = "neighbour")
# dapc1 <- dapc(nephrops_med_gl, nephrops_med_gl$pop, n.pca = 2, n.da = 2)
# 
# # Extract percent of genetic variance is explained by each axis
# (percent <- dapc1$eig/sum(dapc1$eig)*100)
# 
# # Plot DAPC
# dapc1_plt <- plot_pca(
#   dapc1$ind.coord, percent, nephrops_med_gl, rad_samples_meta, by = "site",
#   cols = cols, axis_lab = "LD", title = "Discriminant analysis of principal components",
#   labsize = 3
# )+
#   # Add additional text to the ggplot
#   annotate("text", x = 4.5, y = -2, label = "n.pca = 3", size = 3)
# dapc1_plt
# # plot_pca(dapc1$ind.coord, percent, nephrops_med_gl, rad_samples_meta, by = "sex", axis_lab = "LD")

# # Assess number of clusters using AICc
# AICc <- snapclust.choose.k(gl2gi(nephrops_med_gl), max = nPop(nephrops_med_gl), IC = "AICc")
# 
# # Plot AICc
# AICc_plt <- ggplot(data = data.frame(AICc = AICc, K = 1:length(AICc)))+
#   geom_point(aes(x=K, y=AICc), shape = 19, size = 2.5, colour = "black")+
#   geom_line(aes(x=K, y=AICc))+
#   scale_x_continuous(breaks = 1:length(AICc))+
#   ylab("AICc")+
#   xlab("K")+
#   # ggtitle("Assess number of clusters using AICc")+
#   theme(
#     panel.grid.minor.x = element_blank(),
#     axis.title.x = element_text(face = "italic"),
#     plot.title = element_text(hjust=0.5, size=15)
#   )
# AICc_plt
# 
# # K
# K <- 2
# 
# # Impute missing data then run find.clusters()
# grp <- nephrops_med_gl |>
#   gl.impute(x = _, method = "neighbour") |>
#   find.clusters(x = _, n.clust = K, n.pca = 200)
# 
# # Run snapclust
# snapclust_res <- nephrops_med_gl |>
#   gl.impute(x = _, method = "neighbour") |>
#   gl2gi(x = _) |>
#   snapclust(x = _, k = K, pop.ini = grp$grp)
# compoplot(snapclust_res, show.lab = TRUE, cex.names = 0.5, legend = FALSE)

#--------------#
# IBD with Med Sites Only ####
#--------------#

# Read in sites coordinates
sites <- read_csv("./data/site_coordinates.csv") |>
  arrange(.data = _, Code) |>
  filter(.data = _, Code != "Cly")
sites$Code == popNames(nephrops_med_gl)

# Geographic distances
library(marmap)
bathy <- raster::raster("data/MS_bathy_5m_lonlat.tif") |> raster::crop(x = _, y = raster::extent(0,28,35,50))
trans1 <- trans.mat(as.bathy(bathy), min.depth = -10, max.depth = NULL)
# lc_paths <- lc.dist(trans1, subset(sites, select = c("Lon","Lat")), res = "path")
# plot(bathy); lapply(lc_paths, lines, col = "orange", lwd = 2, lty = 1)
lc_dist <- lc.dist(trans1, subset(sites, select = c("Lon","Lat")), res = "dist")

# Genetic distances
Fst_med <- as.dist(as.matrix(Fst)[-c(7),-c(7)])
plot(lc_dist, Fst_med)

# Isolation-by-distance Mantel test
library(ade4)
IBD = mantel.rtest(lc_dist, Fst_med, nrepet = 999)
IBD

# Plot geographic versus genetic distances
Fst_med_df <- data.frame(
  Fst = as.vector(Fst_med),
  Geo = as.vector(lc_dist)
)
IBD_plt <- ggplot(data = Fst_med_df)+
  geom_point(aes(x = Geo, y = Fst), size = 3.5, shape = 21, colour = "white", fill = "black", stroke = 0.1)+
  scale_y_continuous(limits = c(0, 0.030))+
  annotate("text", x = 200, y = 0.029, size = 3, label = expression(italic("r")^2 ~ "= 0.32," ~ italic("p") ~ "= 0.23" ))+
  xlab("Geographic distances (km)")+
  ylab(fst.label)+
  ggtitle("Isolation-by-distance test: Mediterranean sites")+
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.2),
    panel.background = element_rect(fill = "#f0f0f0"),
  )
IBD_plt

#--------------#
# Figure 3 ####
#--------------#

# Plot layout
plt_fig3 = list(
  pca2_plt+ labs(tag = "A")+ theme(
    axis.title = element_text(size=10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size=12),
  ),
  IBD_plt+ labs(tag = "B")
)
figure3 <- wrap_plots(plt_fig3)

# Export
ggsave(plot = figure3, filename = "figure3.pdf", width = 11, height = 5)
