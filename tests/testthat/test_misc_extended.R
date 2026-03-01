context("misc extended")

# ---- parse_region ----
test_that("parse_region parses valid region string", {
  result <- parse_region("chr1:100-200")
  expect_s3_class(result, "data.frame")
  expect_equal(result$chrom, "1")
  expect_equal(result$start, 100L)
  expect_equal(result$end, 200L)
})

test_that("parse_region handles X chromosome", {
  result <- parse_region("chrX:500-1000")
  expect_equal(result$chrom, "X")
  expect_equal(result$start, 500L)
  expect_equal(result$end, 1000L)
})

test_that("parse_region errors on invalid format", {
  expect_error(parse_region("1:100-200"), "format must be")
  expect_error(parse_region("chr1-100-200"), "format must be")
  expect_error(parse_region("chr1:abc-200"), "format must be")
})

test_that("parse_region returns non-string input unchanged", {
  df <- data.frame(chrom = 1, start = 100, end = 200)
  result <- parse_region(df)
  expect_identical(result, df)
})

test_that("parse_region returns non-single-string input unchanged", {
  input <- c("chr1:100-200", "chr2:300-400")
  result <- parse_region(input)
  expect_identical(result, input)
})

# ---- parse_variant_id ----
test_that("parse_variant_id parses single variant with chr prefix", {
  result <- parse_variant_id("chr1:12345:A:G")
  expect_equal(result$chrom, 1L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
})

test_that("parse_variant_id parses single variant without chr prefix", {
  result <- parse_variant_id("5:12345:A:G")
  expect_equal(result$chrom, 5L)
  expect_equal(result$pos, 12345L)
  expect_equal(result$A2, "A")
  expect_equal(result$A1, "G")
  conv <- attr(result, "convention")
  expect_false(conv$has_chr)
})

test_that("parse_variant_id parses multiple variants", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T", "chr3:300:G:A")
  result <- parse_variant_id(ids)
  expect_equal(nrow(result), 3)
  expect_equal(result$chrom, c(1L, 2L, 3L))
  expect_equal(result$pos, c(100L, 200L, 300L))
  expect_equal(result$A2, c("A", "C", "G"))
  expect_equal(result$A1, c("G", "T", "A"))
})

# ---- detect_variant_convention ----
test_that("detect_variant_convention detects chr prefix and allele separators", {
  conv <- pecotmr:::detect_variant_convention(c("chr1:100:A:G", "chr2:200:C:T"))
  expect_true(conv$has_chr)
  expect_equal(conv$allele_sep, ":")
  expect_false(conv$has_build)

  conv2 <- pecotmr:::detect_variant_convention(c("1_100_A_G", "2_200_C_T"))
  expect_false(conv2$has_chr)
  expect_equal(conv2$allele_sep, "_")

  conv3 <- pecotmr:::detect_variant_convention(c("chr1:100:A:G:b38"))
  expect_true(conv3$has_build)

  # Mixed format: colon between chrom:pos, underscore between alleles
  conv4 <- pecotmr:::detect_variant_convention(c("chr1:100_A_G"))
  expect_true(conv4$has_chr)
  expect_equal(conv4$allele_sep, "_")

  # Mixed format without chr prefix
  conv5 <- pecotmr:::detect_variant_convention(c("1:100_A_G"))
  expect_false(conv5$has_chr)
  expect_equal(conv5$allele_sep, "_")
})

# ---- normalize_variant_id ----
test_that("normalize_variant_id normalizes various formats", {
  # Underscore to canonical
  expect_equal(normalize_variant_id("1_100_A_G"), "chr1:100:A:G")
  # Already canonical
  expect_equal(normalize_variant_id("chr1:100:A:G"), "chr1:100:A:G")
  # Without chr prefix, default adds it
  expect_equal(normalize_variant_id("1:100:A:G"), "chr1:100:A:G")
  # With chr_prefix = FALSE
  expect_equal(normalize_variant_id("chr1:100:A:G", chr_prefix = FALSE), "1:100:A:G")
  # With build suffix
  expect_equal(normalize_variant_id("chr1:100:A:G:b38"), "chr1:100:A:G")
  # Mixed format to canonical
  expect_equal(normalize_variant_id("chr1:100_A_G"), "chr1:100:A:G")
  # Convention-driven output preserves user format
  conv <- pecotmr:::detect_variant_convention(c("chr1:100_A_G"))
  expect_equal(normalize_variant_id("1:200:C:T", convention = conv), "chr1:200_C_T")
})

# ---- format_variant_id convention-driven ----
test_that("format_variant_id uses convention parameter automatically", {
  # Canonical convention (default)
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G"), "chr1:100:A:G")

  # Convention-driven: mixed format
  conv_mixed <- list(has_chr = TRUE, allele_sep = "_")
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", convention = conv_mixed), "chr1:100_A_G")

  # Convention-driven: no chr, underscore alleles
  conv_nochr <- list(has_chr = FALSE, allele_sep = "_")
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", convention = conv_nochr), "1:100_A_G")

  # Convention overrides explicit chr_prefix parameter
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", chr_prefix = FALSE, convention = conv_mixed), "chr1:100_A_G")
})

# ---- variant_id_to_df ----
test_that("variant_id_to_df handles colon-separated format", {
  ids <- c("1:100:A:G", "2:200:C:T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
  expect_equal(result$A2, c("A", "C"))
  expect_equal(result$A1, c("G", "T"))
})

test_that("variant_id_to_df handles underscore-separated format", {
  ids <- c("1:100_A_G", "2:200_C_T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$A2, c("A", "C"))
})

test_that("variant_id_to_df strips chr prefix", {
  ids <- c("chr1:100:A:G", "chr2:200:C:T")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(result$chrom, c(1L, 2L))
})

test_that("variant_id_to_df handles data.frame input with named columns", {
  df <- data.frame(chrom = c("chr1", "2"), pos = c(100, 200),
                   A2 = c("A", "C"), A1 = c("G", "T"))
  suppressWarnings(result <- pecotmr:::variant_id_to_df(df))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$pos, c(100L, 200L))
})

test_that("variant_id_to_df handles 5-part IDs with build suffix", {
  ids <- c("chr1:100:A:G:b38", "chr2:200:T:C")
  result <- pecotmr:::variant_id_to_df(ids)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("chrom", "pos", "A2", "A1"))
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$A2, c("A", "T"))
  expect_equal(result$A1, c("G", "C"))
})

test_that("variant_id_to_df handles mixed 4/5-part IDs", {
  ids <- c("1:100:A:G", "chr2:200:T:C:b38", "3:300:G:A:b37")
  suppressWarnings(result <- pecotmr:::variant_id_to_df(ids))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)
  expect_equal(result$A1, c("G", "C", "A"))
})

# ---- get_nested_element ----
test_that("get_nested_element retrieves deeply nested values", {
  nested <- list(a = list(b = list(c = 42)))
  expect_equal(get_nested_element(nested, c("a", "b", "c")), 42)
})

test_that("get_nested_element returns NULL for NULL name_vector", {
  nested <- list(a = 1)
  expect_null(get_nested_element(nested, NULL))
})

test_that("get_nested_element errors on missing element", {
  nested <- list(a = list(b = 1))
  expect_error(get_nested_element(nested, c("a", "x")), "Element not found")
})

test_that("get_nested_element handles single level", {
  nested <- list(a = "hello")
  expect_equal(get_nested_element(nested, "a"), "hello")
})

test_that("get_nested_element skips empty strings", {
  nested <- list(a = list(b = 99))
  expect_equal(get_nested_element(nested, c("", "a", "b")), 99)
})

# ---- region_to_df ----
test_that("region_to_df converts underscore-separated region IDs", {
  ids <- c("1_100_200", "2_300_400")
  result <- region_to_df(ids)
  expect_equal(nrow(result), 2)
  expect_equal(result$chrom, c(1L, 2L))
  expect_equal(result$start, c(100L, 300L))
  expect_equal(result$end, c(200L, 400L))
})

test_that("region_to_df handles chr prefix", {
  ids <- c("chr1_100_200")
  result <- region_to_df(ids)
  expect_equal(result$chrom, 1L)
})

test_that("region_to_df allows custom column names", {
  ids <- c("1_100_200")
  result <- region_to_df(ids, colnames = c("chr", "begin", "finish"))
  expect_true(all(c("chr", "begin", "finish") %in% colnames(result)))
})

# ---- z_to_pvalue ----
test_that("z_to_pvalue returns correct p-values", {
  expect_equal(z_to_pvalue(0), 1.0, tolerance = 1e-10)
  expect_true(z_to_pvalue(1.96) < 0.05)
  expect_true(z_to_pvalue(-1.96) < 0.05)
  expect_equal(z_to_pvalue(1.96), z_to_pvalue(-1.96), tolerance = 1e-10)
})

test_that("z_to_pvalue handles vector input", {
  z <- c(0, 1, 2, 3)
  p <- z_to_pvalue(z)
  expect_length(p, 4)
  expect_true(all(diff(p) < 0))  # p-values should decrease as z increases
})

test_that("z_to_pvalue handles extreme values", {
  expect_true(z_to_pvalue(10) < 1e-20)
  expect_true(z_to_pvalue(40) >= 0)  # Should not be negative
})

# ---- wald_test_pval ----
test_that("wald_test_pval computes correct p-values", {
  # A large t-value should give small p-value
  pval <- wald_test_pval(beta = 5, se = 1, n = 100)
  expect_true(pval < 0.001)

  # Zero beta should give p-value of 1
  pval_zero <- wald_test_pval(beta = 0, se = 1, n = 100)
  expect_equal(pval_zero, 1.0, tolerance = 1e-10)
})

test_that("wald_test_pval handles vector inputs", {
  betas <- c(0, 1, 2, 5)
  ses <- c(1, 1, 1, 1)
  pvals <- wald_test_pval(betas, ses, n = 100)
  expect_length(pvals, 4)
  expect_true(pvals[1] > pvals[4])
})

# ---- compute_LD ----
test_that("compute_LD sample method produces valid correlation matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "sample")
  expect_equal(nrow(R), 4)
  expect_equal(ncol(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
  expect_true(all(R >= -1 & R <= 1))
})

test_that("compute_LD population method produces valid matrix", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "population")
  expect_equal(nrow(R), 4)
  expect_equal(unname(diag(R)), rep(1, 4))
  expect_true(isSymmetric(R))
})

test_that("compute_LD with NA genotypes and sample method", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  X[1, 1] <- NA
  X[5, 3] <- NA
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "sample")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 4))
})

test_that("compute_LD with NA genotypes and population method", {
  set.seed(42)
  X <- matrix(sample(0:2, 200, replace = TRUE), nrow = 50)
  X[1, 1] <- NA
  X[5, 3] <- NA
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "population")
  expect_true(all(!is.na(R)))
  expect_equal(unname(diag(R)), rep(1, 4))
})

test_that("compute_LD errors on NULL input", {
  expect_error(compute_LD(NULL), "X must be provided")
})

test_that("compute_LD with trim_samples and population method", {
  set.seed(42)
  # 13 samples, not multiple of 4
  X <- matrix(sample(0:2, 52, replace = TRUE), nrow = 13)
  colnames(X) <- paste0("rs", 1:4)

  R <- compute_LD(X, method = "population", trim_samples = TRUE)
  expect_equal(nrow(R), 4)
})

test_that("compute_LD with monomorphic SNP", {
  X <- matrix(c(rep(1, 50), sample(0:2, 50, replace = TRUE)), nrow = 50)
  colnames(X) <- c("mono", "poly")

  R <- compute_LD(X, method = "sample")
  expect_equal(unname(diag(R)), c(1, 1))
  # Correlation with monomorphic should be 0 (set by function)
  expect_true(R[1, 2] == 0 || abs(R[1, 2]) < 1e-10)
})

test_that("compute_LD sample vs population differ but are close", {
  set.seed(42)
  X <- matrix(sample(0:2, 500, replace = TRUE), nrow = 100)
  colnames(X) <- paste0("rs", 1:5)

  R_sample <- compute_LD(X, method = "sample")
  R_pop <- compute_LD(X, method = "population")

  # Should be different but close
  expect_false(identical(R_sample, R_pop))
  expect_true(max(abs(R_sample - R_pop)) < 0.1)
})

# ---- filter_X_with_Y ----
test_that("filter_X_with_Y drops variants monomorphic due to Y missingness", {
  X <- matrix(c(0, 0, 1, 2, 0, 1, 1, 2), nrow = 4, ncol = 2)
  Y <- matrix(c(1, NA, NA, 1, 1, 1, 1, NA), nrow = 4, ncol = 2)
  colnames(Y) <- c("ctx1", "ctx2")
  rownames(X) <- rownames(Y) <- paste0("s", 1:4)

  # X column 1 becomes c(0, 0) for ctx1 NAs (s2, s3) -> not monomorphic, still c(0, 2)
  # but after filter_X:
  result <- pecotmr:::filter_X_with_Y(X, Y, missing_rate_thresh = 1, maf_thresh = 0)
  expect_true(is.matrix(result))
})

# ---- filter_Y ----
test_that("filter_Y removes columns with too few non-missing values", {
  Y <- matrix(c(1, NA, NA, NA, 1, 2, 3, 4), nrow = 4)
  result <- pecotmr:::filter_Y(Y, n_nonmiss = 3)
  # Column 1 has only 1 non-missing, column 2 has 4 - column 1 dropped
  # When only 1 column left, may become vector; check length
  expect_true(length(result$Y) >= 3)
})

test_that("filter_Y removes all-NA rows from matrix", {
  Y <- matrix(c(NA, NA, 1, 2, NA, NA, 3, 4), nrow = 4)
  result <- pecotmr:::filter_Y(Y, n_nonmiss = 1)
  expect_true(nrow(result$Y) < 4)
})

# ---- z_to_beta_se (internal) ----
test_that("z_to_beta_se produces correct conversions", {
  z <- c(2.0, -1.0)
  maf <- c(0.3, 0.1)
  n <- 10000

  result <- pecotmr:::z_to_beta_se(z, maf, n)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("beta", "se", "maf") %in% names(result)))
  expect_equal(nrow(result), 2)
  # Sign of beta should match sign of z
  expect_true(result$beta[1] > 0)
  expect_true(result$beta[2] < 0)
})

test_that("z_to_beta_se errors on mismatched lengths", {
  expect_error(pecotmr:::z_to_beta_se(c(1, 2), c(0.3), 1000), "same length")
})

test_that("z_to_beta_se adjusts MAF > 0.5", {
  result <- pecotmr:::z_to_beta_se(c(1.0), c(0.7), 1000)
  expect_equal(result$maf, 0.3)
})

# ---- pval_cauchy (internal) ----
test_that("pval_cauchy combines p-values", {
  pvals <- c(0.01, 0.05, 0.5)
  combined <- pecotmr:::pval_cauchy(pvals)
  expect_true(combined > 0 && combined < 1)
})

test_that("pval_cauchy handles NAs with na.rm", {
  pvals <- c(0.01, NA, 0.05)
  combined <- pecotmr:::pval_cauchy(pvals, na.rm = TRUE)
  expect_true(!is.na(combined))
})

test_that("pval_cauchy with very small p-values", {
  pvals <- c(1e-20, 1e-15)
  combined <- pecotmr:::pval_cauchy(pvals)
  expect_true(combined < 1e-10)
})

# ---- pval_acat (internal) ----
test_that("pval_acat returns single p-value unchanged", {
  expect_equal(pecotmr:::pval_acat(0.05), 0.05)
})

test_that("pval_acat combines multiple p-values", {
  pvals <- c(0.001, 0.01, 0.1)
  combined <- pecotmr:::pval_acat(pvals)
  expect_true(combined > 0 && combined < 1)
})

# ---- matxMax (internal) ----
test_that("matxMax finds location of maximum", {
  mtx <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  result <- pecotmr:::matxMax(mtx)
  expect_equal(result[1], 2)  # row
  expect_equal(result[2], 3)  # col
})

# ---- lbf_to_alpha ----
test_that("lbf_to_alpha converts matrix correctly", {
  lbf <- matrix(c(-0.5, 1.2, 0.3, 0.7, -1.1, 0.4), nrow = 2)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  # Each row should sum to ~1
  expect_equal(rowSums(result), c(1, 1), tolerance = 1e-10)
})

test_that("lbf_to_alpha handles single column", {
  lbf <- matrix(c(0.5, 1.0), ncol = 1)
  colnames(lbf) <- "v1"
  result <- lbf_to_alpha(lbf)
  expect_equal(ncol(result), 1)
})

test_that("lbf_to_alpha handles all-zero row", {
  lbf <- matrix(c(0, 0, 0), nrow = 1)
  colnames(lbf) <- paste0("v", 1:3)
  result <- lbf_to_alpha(lbf)
  expect_true(all(result == 0))
})

# ---- find_data ----
test_that("find_data retrieves from nested list at depth 2", {
  x <- list(a = list(val = 42), b = list(val = 99))
  result <- find_data(x, c(2, "val"))
  expect_equal(result, c(42, 99))
})

test_that("find_data returns list at depth 0", {
  x <- list(a = 1, b = 2)
  result <- find_data(x, c(0))
  expect_type(result, "list")
})

# ---- format_variant_id ----
test_that("format_variant_id constructs canonical IDs", {
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G"), "chr1:100:A:G")
  expect_equal(pecotmr:::format_variant_id(c(1, 2), c(100, 200), c("A", "C"), c("G", "T")),
               c("chr1:100:A:G", "chr2:200:C:T"))
  # Without chr prefix
  expect_equal(pecotmr:::format_variant_id(1, 100, "A", "G", chr_prefix = FALSE), "1:100:A:G")
  # Handles chr-prefixed input chrom
  expect_equal(pecotmr:::format_variant_id("chr1", 100, "A", "G"), "chr1:100:A:G")
})

# ---- is_zero_variance (internal) ----
test_that("is_zero_variance detects constant vector", {
  expect_true(pecotmr:::is_zero_variance(c(5, 5, 5, 5)))
  expect_false(pecotmr:::is_zero_variance(c(1, 2, 3)))
})

# ---- mean_impute (internal) ----
test_that("mean_impute fills NAs with column means", {
  X <- matrix(c(1, 2, NA, 4, NA, 6), nrow = 3)
  result <- pecotmr:::mean_impute(X)
  expect_true(all(!is.na(result)))
  expect_equal(result[3, 1], 1.5)  # mean of 1, 2
  expect_equal(result[2, 2], 5)    # mean of 4, 6
})
