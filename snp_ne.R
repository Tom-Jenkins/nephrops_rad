# ============================ #
#
# Nephrops RAD Analysis R Script 2023
#
# SNP Effective Population Size
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
library(stringr)
library(dartR)
library(qvalue)
library(OutFLANK)

# Import QC SNPs
load("./data/nephrops_snps.RData"); nephrops_gl

# ----------------- #
# Ne
# ----------------- #

# Create an Adriatic versus Pomo Pit data set
nephrops_ne_gl <- nephrops_gl |>
  gl.drop.pop(x = _, pop.list = c("9I","Aeg","Cly"), mono.rm = TRUE)
nephrops_ne_gl$pop <- nephrops_ne_gl$pop |>
  str_replace_all(string = _, "17I|18II|Anc|Cgg", "Adriatic") |>
  str_replace_all(string = _, "Pom1|Pom2|Pom3", "Pomo") |>
  as.factor(x = _)
nephrops_ne_gl; summary(nephrops_ne_gl$pop)

# Check for monomorphic loci
gl.report.monomorphs(nephrops_ne_gl)

# Search for outlier loci using OutFLANK
outflank_results <- gl.outflank(nephrops_ne_gl)$outflank$results

# Number of outliers (no outliers detected)
length(which(outflank_results$OutlierFlag == TRUE))

# Run NeEstimator V2.1
ne <- gl.LDNe(nephrops_ne_gl,
              outfile = "nephropsLDNE.txt",
              outpath = ".",
              neest.path = "./NeEstimatorv2x1/",
              critical = c(0.01, 0.02, 0.05),
              singleton.rm = TRUE,
              mating = "random",
)

# Show results
ne

# ----------------- #
# NeEstimator Manual Notes
# ----------------- #

# Confidence Intervals:
# A new method (Jones et al, 2016) is implemented for calculating
# confidence intervals by jackknifing over individuals rather than loci, as
# originally suggested by Waples and Do (2008). This new method
# accounts for pseudoreplication due to physical linkage and overlapping
# pairs of loci being compared. Parametric CIs are still reported also, but
# they are much too narrow with large numbers of loci, so the new
# jackknifed CIs should be used with all large SNP datasets.

# The user is not likely to see large differences between the new jackknife CIs
# and parametric CIs for small numbers of loci, but we recommend the general 
# use of the jackknife CIs, particularly when the number of loci is large (>100).
# For publication, users should describe which types of confidence intervals 
# (i.e. parametric, or jackknifed) are reported alongside Ne estimates.

# Confidence intervals for all methods assume that the loci assort independently.
# As noted above, if some pairs of loci are physically linked, the data contain less
# information than assumed and the confidence intervals will be too narrow.
# Physically linked loci would in general also downwardly bias estimates of Ne
# from the LD method. The consequences of these departures from the
# assumption of independent assortment have rigorously evaluated by Waples et
# al (2016).

# Rare Alleles:
# A new option is provided for screening out rare alleles. Current options
# involve choosing one or more fixed PCrit values that specify maximum
# allowable allele frequencies. The problem is that the consequences of
# any fixed PCrit depends on the sample size, which can vary among
# populations in an input file and, within populations, among pairs of loci
# because of varying degrees of missing data. The latter is expected to be
# an important issue in large SNP datasets obtained using NGS methods.
# The user can now select a new option, which only removes singleton
# alleles—those that occur in only a single copy in one heterozygote. These
# singleton alleles contribute the most to upward bias in Ne^ for the LD
# method. This new option in effect allows PCrit to vary across locus pairs.

# Negative / infinite values interpretation:
# Because the actual contribution of sampling error is a random variable, it can
# be smaller than the expected value, and when that happens subtracting the
# expected contribution can produce a negative estimate of adjusted Fˆ or rˆ2 ,
# which in turn produces a negative estimate of Ne. This also can occur with
# unbiased estimators of FST or genetic distance.

# The usual interpretation in this case is that the estimate of Ne is infinity – that
# is, there is no evidence for variation in the genetic characteristic caused by
# genetic drift due to a finite number of parents — it can all be explained by
# sampling error.

# In the confidence intervals (CIs), such values are reported as ‘Infinite,’ meaning
# that the confidence interval includes infinity. However, the point estimates of Ne
# are reported even if they are negative (in accessory output files), as in some
# applications this information can be useful.

# For example, say you have several replicate samples from the same population
# and use each sample to estimate Ne. An overall estimate of Ne can be obtained
# by taking the harmonic mean of the separate estimates, even if they are
# negative. You will get an approximately, but not exactly, correct answer if you
# replace negative estimates of Ne with infinity before taking a harmonic mean.
# This issue is discussed in Waples and Do (2010).

