# ============================ #
#
# Nephrops RAD Analysis R Script 2023
#
# SNP QC and Filtering
#
# Data source:
# # ./data/populations.snps.vcf
# VCF file output from Stacks populations v2.53.
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(vcfR)
library(adegenet)
library(poppr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(dartR)
# dartR::gl.install.vanilla.dartR() # install dartR dependencies
library(snpStats) # BiocManager::install("snpStats")
library(qvalue) # BiocManager::install("qvalue")
library(OutFLANK) # devtools::install_github("whitlock/OutFLANK")

# ----------------- #
# Data Import
# ----------------- #

# Import vcf file
vcf <- read.vcfR("./data/populations.snps.vcf", verbose = FALSE)
vcf

# Check for biallelic SNPs
summary(is.biallelic(vcf))

# Only keep biallelic SNPs
vcf <- vcf[is.biallelic(vcf) == TRUE]
vcf

# Convert vcf to genlight object
nephrops <- vcfR2genind(vcf)
nephrops

# Print individual names
indNames(nephrops)

# Add population labels
popID <- read.delim("./data/popmap_filt.txt", header = FALSE)
nrow(popID) == length(indNames(nephrops))
nephrops$pop <- as.factor(popID$V2)
summary(nephrops$pop)

# ----------------- #
# Individual Missing Data QC
# ----------------- #

# Missing data for individuals
indmiss <- propTyped(nephrops, by = "ind")

# Plot using barplot
barplot(indmiss, ylim = c(0,1), ylab = "Complete genotypes (proportion)", las = 2)
indmiss[ which(indmiss < 0.70) ] # print ind with less than 70% complete genotypes

# Remove individuals with >30% missing data
nephrops <- missingno(nephrops, type = "geno", cutoff = 0.30)
summary(nephrops$pop)

# Missing data for loci
locmiss <- propTyped(nephrops, by = "ind")
summary(locmiss)

# Remove SNPs with >30% missing data
nephrops <- missingno(nephrops, type = "loci", cutoff = 0.30)

# Filter loci that depart from HWE
nephrops_gl <- gi2gl(nephrops)
nephrops_gl <- gl.filter.hwe(nephrops_gl, n.pop.threshold = 3, min_sample_size = 9)

# Calculate and filter loci that are in LD
ld_by_pop <- gl.report.ld.map(nephrops_gl, ind.limit = 9, maf = 3) # maf > 1 == mac
nephrops_gl <- gl.filter.ld(nephrops_gl, ld_by_pop, threshold = 0.5, pop.limit = 3)

# Check for loci that are all NA in a single site
gl.report.heterozygosity(nephrops_gl, plot.out = FALSE)$all_NALoc

# Remove any loci that are all NA in one site
nephrops_gl <- gl.filter.allna(nephrops_gl, by.pop = TRUE)

# Check for monomorphic loci
gl.report.monomorphs(nephrops_gl)

# Check minor allele frequency
gl.report.maf(nephrops_gl)

# Keep loci with a minor allele count >= 5
nephrops_gl <- gl.filter.maf(nephrops_gl, threshold = 5)

# ----------------- #
# Neutral Loci Check
# ----------------- #

# Perform outFlank
outflank_results <- gl.outflank(nephrops_gl,
                                LeftTrimFraction = 0.05,
                                RightTrimFraction = 0.20)
# Number of outliers
(numOutliers <- length(which(outflank_results$OutlierFlag == TRUE)))

# If there are outliers, create a neutral SNP genlight object
if (numOutliers != 0) {
  outliers <- outflank_results[which(outflank_results$OutlierFlag == TRUE),]$LocusName
  nephrops_gl <- gl.drop.loc(nephrops_gl, loc.list = outliers)
}

#--------------#
# Export Data
#--------------#

# Export genlight object as .RData file
nephrops_gl
summary(nephrops_gl$pop)
save(nephrops_gl, file = "./data/nephrops_snps.RData")

# Export genlight object as .geno file
# nephrops_gl_geno <- as.data.frame(nephrops_gl)
# nephrops_gl_geno[is.na(nephrops_gl_geno)] <- 9
# write.geno(nephrops_gl_geno, output.file = "./Data/nephrops_gl.geno")

# Export genlight object as STRUCTURE file
# gl2structure(nephrops_gl,
#              indNames = indNames(nephrops_gl),
#              addcolumns = factor(nephrops_gl$pop, labels = 1:nPop(nephrops_gl)),
#              exportMarkerNames = TRUE,
#              outfile = "nephrops_gl.struct",
#              outpath = "./Data/"
# )
