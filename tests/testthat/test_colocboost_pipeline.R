context("colocboost_pipeline: qc_method default and pip_cutoff validation")

# ---- qc_method match.arg ----
test_that("qc_regional_data resolves qc_method default to 'dentist'", {
  # With no data, it should return the region_data unchanged
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  result <- pecotmr:::qc_regional_data(region_data)
  expect_type(result, "list")
})

test_that("qc_regional_data accepts explicit qc_method = 'dentist'", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  result <- pecotmr:::qc_regional_data(region_data, qc_method = "dentist")
  expect_type(result, "list")
})

test_that("qc_regional_data accepts explicit qc_method = 'slalom'", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  result <- pecotmr:::qc_regional_data(region_data, qc_method = "slalom")
  expect_type(result, "list")
})

test_that("qc_regional_data rejects invalid qc_method", {
  region_data <- list(individual_data = NULL, sumstat_data = NULL)
  expect_error(
    pecotmr:::qc_regional_data(region_data, qc_method = "invalid"),
    "arg"
  )
})

# ---- pip_cutoff_to_skip_ind validation ----
test_that("pip_cutoff scalar is recycled for individual contexts", {
  # Create individual_data with 3 real-ish contexts
  set.seed(42)
  n <- 10; p <- 5
  make_ctx <- function() {
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("var", 1:p)
    Y <- matrix(rnorm(n * 2), n, 2)
    colnames(Y) <- paste0("gene", 1:2)
    list(X = X, Y = Y, maf = runif(p, 0.05, 0.5))
  }
  ctx <- make_ctx()
  individual_data <- list(
    residual_Y = list(ctx1 = ctx$Y, ctx2 = ctx$Y, ctx3 = ctx$Y),
    residual_X = list(ctx1 = ctx$X, ctx2 = ctx$X, ctx3 = ctx$X),
    maf = list(ctx1 = ctx$maf, ctx2 = ctx$maf, ctx3 = ctx$maf)
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  # Scalar 0 (no PIP check) should be recycled and run without error
  result <- pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = 0)
  expect_type(result, "list")
})

test_that("pip_cutoff wrong length errors for individual contexts", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  individual_data <- list(
    residual_Y = list(ctx1 = Y, ctx2 = Y, ctx3 = Y),
    residual_X = list(ctx1 = X, ctx2 = X, ctx3 = X),
    maf = list(ctx1 = runif(p), ctx2 = runif(p), ctx3 = runif(p))
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  expect_error(
    pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = c(0, 0)),
    "pip_cutoff_to_skip_ind"
  )
})

test_that("pip_cutoff correct length vector works", {
  set.seed(42)
  n <- 10; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("var", 1:p)
  Y <- matrix(rnorm(n), n, 1)
  colnames(Y) <- "gene1"
  individual_data <- list(
    residual_Y = list(ctx1 = Y, ctx2 = Y),
    residual_X = list(ctx1 = X, ctx2 = X),
    maf = list(ctx1 = runif(p, 0.05, 0.5), ctx2 = runif(p, 0.05, 0.5))
  )
  region_data <- list(individual_data = individual_data, sumstat_data = NULL)

  # Length-2 vector for 2 contexts should work
  result <- pecotmr:::qc_regional_data(region_data, pip_cutoff_to_skip_ind = c(0, 0))
  expect_type(result, "list")
})

# ---- is_valid_sumstat_entry ----
test_that("is_valid_sumstat_entry returns TRUE for valid data", {
  ss <- data.frame(z = c(1.5, -2.3, 0.8), n = c(1000, 1000, 1000),
                   variant = c("chr1:100:A:G", "chr1:200:T:C", "chr1:300:G:A"))
  expect_true(pecotmr:::is_valid_sumstat_entry(ss))
})

test_that("is_valid_sumstat_entry returns FALSE for too few variants", {
  ss <- data.frame(z = 1.5, n = 1000, variant = "chr1:100:A:G")
  expect_false(pecotmr:::is_valid_sumstat_entry(ss, min_variants = 2))
})

test_that("is_valid_sumstat_entry returns FALSE for all-NA z-scores", {
  ss <- data.frame(z = c(NA, NA, NA), n = c(1000, 1000, 1000),
                   variant = c("v1", "v2", "v3"))
  expect_false(pecotmr:::is_valid_sumstat_entry(ss))
})

test_that("is_valid_sumstat_entry returns FALSE for NULL input", {
  expect_false(pecotmr:::is_valid_sumstat_entry(NULL))
})

test_that("is_valid_sumstat_entry returns FALSE for zero sample size", {
  ss <- data.frame(z = c(1.5, -2.3), n = c(0, 0), variant = c("v1", "v2"))
  expect_false(pecotmr:::is_valid_sumstat_entry(ss))
})

# ---- filter_valid_sumstats ----
test_that("filter_valid_sumstats removes invalid entries", {
  good <- data.frame(z = c(1.0, 2.0, 3.0), n = c(100, 100, 100),
                     variant = c("v1", "v2", "v3"))
  bad <- data.frame(z = c(NA, NA), n = c(0, 0), variant = c("v1", "v2"))
  sumstats <- list(study_a = good, study_b = bad)
  LD_mat <- list(ld1 = matrix(1, 3, 3))
  LD_match <- c("ld1", "ld1")
  result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match)
  expect_equal(length(result$sumstats), 1)
  expect_equal(names(result$sumstats), "study_a")
})

test_that("filter_valid_sumstats returns NULL when all invalid", {
  bad1 <- data.frame(z = NA, n = 0, variant = "v1")
  bad2 <- data.frame(z = NA, n = 0, variant = "v1")
  sumstats <- list(s1 = bad1, s2 = bad2)
  LD_mat <- list(ld1 = matrix(1))
  LD_match <- c("ld1", "ld1")
  expect_message(
    result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match),
    "No valid"
  )
  expect_null(result)
})

test_that("filter_valid_sumstats keeps all when all valid", {
  make_ss <- function() data.frame(z = c(1, 2, 3), n = c(100, 100, 100),
                                   variant = c("v1", "v2", "v3"))
  sumstats <- list(s1 = make_ss(), s2 = make_ss())
  LD_mat <- list(ld1 = matrix(1, 3, 3))
  LD_match <- c("ld1", "ld1")
  result <- pecotmr:::filter_valid_sumstats(sumstats, LD_mat, LD_match)
  expect_equal(length(result$sumstats), 2)
})
