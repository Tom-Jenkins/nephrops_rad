# ============================ #
#
# Nephrops RAD Analysis R Script 2023
#
# SNP Diversity Stats
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
library(dartR)
library(snpStats)
library(dplyr)
library(tidyr)

# Import QC SNPs
load("./data/nephrops_snps.RData")

# ----------------- #
# Heterozygosity
# ----------------- #

# Plot histogram of loci heterozygosity
hist(dartR::gl.Ho(nephrops_gl))
hist(dartR::gl.He(nephrops_gl))

# Calculate heterozygosity stats per site
het_stats <- gl.report.heterozygosity(nephrops_gl)

# Observed, unbiased expected heterozygosity and Fis per site
het_stats |> select(.data = _, Ho, uHe, FIS) |> round(x = _, digits = 3)

# Number of polymorphic / monomorphic loci per site
het_stats |> select(.data = _, polyLoc, monoLoc)

# Report number of private alleles per site
pa <- gl.report.pa(nephrops_gl, method = "one2rest", loc_names = TRUE, plot.out = FALSE)

# Print number of private and fixed SNPs
pa$table |> select(.data = _, pop1, priv1, fixed)
# private_alleles(gl2gi(nephrops_gl), level = "population", report = "table") |>
#   apply(X = _, MARGIN = 1, FUN = function(row) sum(row >= 1))

# ----------------- #
# Visualise AF of private alleles in Clyde Sea
# ----------------- #

# Extract the ten loci that are private to the Clyde Sea
(clyde_snps <- pa$names_loci$Cly$pop1_pop2_pa)

# Subset SNPs from genlight object
clyde_snps_gl <- nephrops_gl[, clyde_snps]

# Calculate allele frequencies per site
af <- lapply(1:nPop(clyde_snps_gl), function(x) gl.alf(seppop(clyde_snps_gl)[[x]]))
names(af) <- popNames(clyde_snps_gl)

# Extract frequencies for first allele and convert to data.frame
af1 <- lapply(seq_along(af), function(x) af[[x]]$alf1)
names(af1) <- names(af)
af1_df <- list2DF(af1)

# Add loci names
af1_df$loci <- clyde_snps

# Convert data.frame to long format
af1_df <- pivot_longer(af1_df, cols = !loci, names_to = "site", values_to = "freq")

# Plot barplot
ggplot(data=af1_df, aes(x=site, y=freq, fill=site))+
  geom_bar(stat="identity", colour="black", size=0.5)+
  facet_wrap(~loci, scales="free")+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))+
  ylab("Allele frequency")+
  xlab("Sampling site")+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black", size=6),
    axis.ticks.x = element_blank(),
    axis.title = element_text(colour="black", size=15),
    strip.text = element_text(colour="black", size=14), 
    panel.background = element_rect(fill="white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth=0.5),
    plot.title = element_text(hjust = 0.5, size=18), # title centered
    legend.title = element_blank(),
    legend.text = element_text(size=15),
    legend.position = "top",
    legend.justification = "centre"
  )
