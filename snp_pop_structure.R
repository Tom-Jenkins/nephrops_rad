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
load("./Data/nephrops_snps.RData"); nephrops_gl

# Import RAD samples metadata
rad_samples_meta <- read.csv("./Data/growth_data_rad_samples_only.csv")

# Assess sex ratio per site
rad_samples_meta |>
  filter(.data = _, Ind_ID %in% indNames(nephrops_gl)) |>
  count(x = _, Site, Sex)

#--------------#
# FST
#--------------#

# Compute pairwise Fst
Fst <- genet.dist(gl2gi(nephrops_gl), method = "WC84")

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

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] <- 0

# Round digits
fst.df$Fst <- round(fst.df$Fst, digits = 2)

# Fst italic label
fst.label <- expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid <- max(fst.df$Fst) / 2

# Plot heatmap
fst_heat <- ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="white", size = 2.5)+
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
# Perform PCA with Clyde Sea
#--------------#

# Impute missing data then run PCA
pca1 <- nephrops_gl |>
  gl.impute(x = _, method = "neighbour") |>
  # tab(x = _, NA.method = "mean") |>
  glPca(x = _, scale = TRUE, center = TRUE, nf = 3)

# Extract percent of genetic variance is explained by each axis
pca1$eig/sum(pca1$eig)*100

# PCA colours
cols = c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#d6604d"
         ,"#2171b5","#2171b5","#2171b5")

# Plot PCA
(pca_plt <- plot_pca(pca1, nephrops_gl, rad_samples_meta, by = "site", cols = cols))
# pca_plt <- plot_pca(pca1, nephrops_gl, rad_samples_meta, by = "sex", axes = c(1,2))

#--------------#
# Perform Snapclust with Clyde Sea
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
# write.table(mem_prob, file = "./Data/snapclust_membership_prob.txt", row.names = FALSE)

#--------------#
# Figure 2
#--------------#

# Read in mapmixture coordinates file
coords <- read.csv("Data/mapmixture_coordinates.csv")

# Run mapmixture
library(mapmixture)
mapmixture_plt <- mapmixture(
  admixture_df = mem_prob,
  coords_df = coords,
  cluster_names = c("Cluster 1", "Cluster 2"),
  cluster_cols = c("#2166AC", "#D6604D"),
  boundary = c(xmin=-5.5, xmax=26.5, ymin=36.5, ymax=58),
  crs = 3035,
  arrow_position = "bl",
  arrow_size = 1.5,
  scalebar_position = "bl",
  scalebar_size = 1.3,
  axis_text_size = 10,
  axis_title_size = 12
)+
  theme(
    legend.key.size = unit(0.7, "cm"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    legend.spacing.y = unit(0.2, "cm")
  )+
  guides(fill = guide_legend(byrow = TRUE))

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
# Perform PCA with Med Sites Only
#--------------#

# # Data set for PCA (remove Clyde Sea and subsequent monomorphic SNPs)
# nephrops_med_gl <- nephrops_gl |>
#   gl.drop.pop(x = _, pop.list = c("Cly"), mono.rm = TRUE)
# 
# # Impute missing data then run PCA
# pca2 <- nephrops_med_gl |>
#   gl.impute(x = _, method = "neighbour") |>
#   # tab(x = _, NA.method = "mean") |>
#   dudi.pca(df = _, scannf = FALSE, scale = TRUE, center = TRUE, nf = 3)
# 
# # Extract percent of genetic variance is explained by each axis
# pca2$eig/sum(pca2$eig)*100
# 
# # Plot PCA
# plot_pca(pca2, nephrops_med_gl, rad_samples_meta, by = "site")
# plot_pca(pca2, nephrops_med_gl, rad_samples_meta, by = "sex")

#--------------#
# Perform sNMF with Clyde Sea
#--------------#

# # Run SNMF
# snmf("./Data/nephrops_gl.geno",
#      K = 1:4, # number of K ancestral populations to run
#      repetitions = 10, # ten repetitions for each K
#      entropy = TRUE, # calculate cross-entropy
#      alpha = 100,
#      ploidy = 2,
#      project = "new",
#      seed = 1234)
# 
# # Load snmf project
# snmf1 <- load.snmfProject("./Data/nephrops_gl.snmfProject")
# 
# # Plot cross-entropy results to assess optimal number of K
# # Smaller values of cross-entropy usually mean better runs
# # A plateau usually represents the K that best fits the data
# plot_cross_entropy(snmf1)
# 
# # Plot admixture results for K ancestral populations
# plot_admixture(snmf1, nephrops_gl, K = 2)

