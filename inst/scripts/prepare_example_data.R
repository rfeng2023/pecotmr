#!/usr/bin/env Rscript
#
# prepare_example_data.R
#
# This script recovers the original CLU locus data from git history,
# de-identifies it (synthetic sample names, shifted positions, no gene/locus
# identifiers), runs SuSiE fine-mapping, and saves four example data objects
# for use in package vignettes.
#
# Prerequisites:
#   1. Recover the original files from git:
#      git show b25ccd7:inst/prototype/CLU_gwas.rds > /tmp/CLU_gwas.rds
#      git show b25ccd7:inst/prototype/pseudo_bulk_CLU.ENSG00000120885.rds > /tmp/pseudo_bulk_CLU.rds
#   2. Have susieR installed
#
# Usage:
#   Rscript inst/scripts/prepare_example_data.R

library(susieR)
set.seed(42)

cat("=== Loading recovered data ===\n")
gwas_orig <- readRDS("/tmp/CLU_gwas.rds")
eqtl_orig <- readRDS("/tmp/pseudo_bulk_CLU.rds")

cat(sprintf("GWAS: %d variants\n", nrow(gwas_orig)))
cat(sprintf("eQTL: %d samples x %d variants\n", nrow(eqtl_orig$X), ncol(eqtl_orig$X)))

# ---------------------------------------------------------------------------
# De-identification parameters
# ---------------------------------------------------------------------------
new_chrom <- "chr22"
pos_offset <- 5000000L  # shift all positions

# ---------------------------------------------------------------------------
# 1. De-identify variant information
# ---------------------------------------------------------------------------
cat("\n=== De-identifying variants ===\n")
new_pos <- gwas_orig$POS + pos_offset
new_variant_ids <- paste0(new_chrom, ":", new_pos, ":", gwas_orig$A1, ":", gwas_orig$A2)

cat(sprintf("Original position range: %d - %d\n", min(gwas_orig$POS), max(gwas_orig$POS)))
cat(sprintf("New position range: %d - %d\n", min(new_pos), max(new_pos)))

# ---------------------------------------------------------------------------
# 2. De-identify sample names
# ---------------------------------------------------------------------------
cat("\n=== De-identifying samples ===\n")
n_samples <- nrow(eqtl_orig$X)
new_sample_names <- sprintf("sample_%03d", seq_len(n_samples))

# ---------------------------------------------------------------------------
# 3. Create gwas_sumstats_example
# ---------------------------------------------------------------------------
cat("\n=== Creating gwas_sumstats_example ===\n")
gwas_sumstats_example <- data.frame(
  variant_id = new_variant_ids,
  chrom = new_chrom,
  pos = new_pos,
  A1 = gwas_orig$A1,
  A2 = gwas_orig$A2,
  beta = gwas_orig$beta,
  se = gwas_orig$standard_error,
  z = gwas_orig$Z,
  stringsAsFactors = FALSE
)
cat(sprintf("  %d variants, %d columns\n", nrow(gwas_sumstats_example), ncol(gwas_sumstats_example)))

# ---------------------------------------------------------------------------
# 4. Create eqtl_region_example
# ---------------------------------------------------------------------------
cat("\n=== Creating eqtl_region_example ===\n")
X_new <- eqtl_orig$X
colnames(X_new) <- gsub("_", ":", gsub("^chr8:", paste0(new_chrom, ":"), colnames(X_new)))
# Now shift positions in column names to match new_variant_ids
# Column names are like chr8:27119788:A:G -> chr22:32119788:A:G
colnames(X_new) <- new_variant_ids
rownames(X_new) <- new_sample_names

y_new <- eqtl_orig$y_res
names(y_new) <- new_sample_names

eqtl_region_example <- list(X = X_new, y_res = y_new)
cat(sprintf("  X: %d x %d, y_res: %d\n",
            nrow(eqtl_region_example$X),
            ncol(eqtl_region_example$X),
            length(eqtl_region_example$y_res)))

# ---------------------------------------------------------------------------
# 5. Create gwas_finemapping_example (SuSiE RSS on GWAS z-scores)
# ---------------------------------------------------------------------------
cat("\n=== Running SuSiE RSS for GWAS fine-mapping ===\n")
R <- cor(eqtl_region_example$X)
z_gwas <- gwas_sumstats_example$z
names(z_gwas) <- new_variant_ids

# Use a realistic GWAS sample size; the exact value is not critical for the
# example but is needed for susie_rss to calibrate effect sizes.
gwas_n <- 400000L
gwas_susie <- susie_rss(
  z = z_gwas,
  R = R,
  n = gwas_n,
  L = 10,
  max_iter = 500,
  estimate_prior_variance = TRUE,
  verbose = FALSE
)
names(gwas_susie$pip) <- new_variant_ids

# Trim to essential fields to reduce file size
gwas_susie_trimmed <- list(
  alpha = gwas_susie$alpha,
  pip = gwas_susie$pip,
  V = gwas_susie$V,
  sets = gwas_susie$sets
)
# Store as list where [[1]] is the SuSiE object (matching xqtl_enrichment_wrapper format)
gwas_finemapping_example <- list(gwas_susie_trimmed)
cat(sprintf("  SuSiE converged, %d credible sets\n", length(gwas_susie$sets$cs)))

# ---------------------------------------------------------------------------
# 6. Create qtl_finemapping_example (SuSiE on eQTL individual-level data)
# ---------------------------------------------------------------------------
cat("\n=== Running SuSiE for QTL fine-mapping ===\n")
qtl_susie <- susie(
  X = eqtl_region_example$X,
  y = eqtl_region_example$y_res,
  L = 10,
  max_iter = 500,
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  verbose = FALSE
)
names(qtl_susie$pip) <- new_variant_ids

# Trim to just the fields needed for enrichment analysis
qtl_susie_trimmed <- list(
  alpha = qtl_susie$alpha,
  pip = qtl_susie$pip,
  V = qtl_susie$V,
  sets = qtl_susie$sets
)

# Store as nested list: [[region]][[context]]$susie_result_trimmed + $variant_names
qtl_finemapping_example <- list(
  region_1 = list(
    context_1 = list(
      susie_result_trimmed = qtl_susie_trimmed,
      variant_names = new_variant_ids
    )
  )
)
cat(sprintf("  SuSiE converged, %d credible sets\n", length(qtl_susie$sets$cs)))

# ---------------------------------------------------------------------------
# 7. Save all data objects
# ---------------------------------------------------------------------------
cat("\n=== Saving data objects ===\n")
data_dir <- "data"

save(gwas_sumstats_example,
     file = file.path(data_dir, "gwas_sumstats_example.rda"),
     compress = "xz")
save(eqtl_region_example,
     file = file.path(data_dir, "eqtl_region_example.rda"),
     compress = "xz")
save(gwas_finemapping_example,
     file = file.path(data_dir, "gwas_finemapping_example.rda"),
     compress = "xz")
save(qtl_finemapping_example,
     file = file.path(data_dir, "qtl_finemapping_example.rda"),
     compress = "xz")

# Print file sizes
for (f in list.files(data_dir, pattern = "example", full.names = TRUE)) {
  cat(sprintf("  %s: %s\n", basename(f), format(file.size(f), big.mark = ",")))
}

cat("\nDone! All example data objects created.\n")
