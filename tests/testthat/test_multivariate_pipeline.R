context("multivariate_pipeline")

# ===========================================================================
# Helpers: small synthetic multivariate data
# ===========================================================================

make_mv_data <- function(n = 30, p = 10, r = 3, seed = 42) {
  set.seed(seed)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  rownames(X) <- paste0("sample_", 1:n)
  colnames(X) <- paste0("chr1:", seq(100, by = 100, length.out = p), ":A:G")
  Y <- matrix(rnorm(n * r), nrow = n, ncol = r)
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("cond_", 1:r)
  maf <- colMeans(X) / 2
  list(X = X, Y = Y, maf = maf)
}

# NOTE: In multivariate_analysis_pipeline, the mvsusieR check (requireNamespace)
# runs BEFORE input validation. So when mvsusieR is not installed, the function
# stops with "please install mvsusieR" before ever reaching the X/Y/maf checks.
# Input validation tests therefore require mvsusieR to be installed.

# ===========================================================================
# Input validation tests (require mvsusieR since validation comes after the
# mvsusieR availability check in the pipeline)
# ===========================================================================

test_that("multivariate_analysis_pipeline errors when X is not a numeric matrix", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = as.data.frame(d$X), Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "X must be a numeric matrix"
  )
})

test_that("multivariate_analysis_pipeline errors when Y is not a numeric matrix", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = as.data.frame(d$Y), maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "Y must be a numeric matrix"
  )
})

test_that("multivariate_analysis_pipeline errors when nrow(X) != nrow(Y)", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  Y_short <- d$Y[1:10, , drop = FALSE]
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = Y_short, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "X and Y must have the same number of rows"
  )
})

# ---------- maf validation -------------------------------------------------

test_that("multivariate_analysis_pipeline errors when maf has wrong length", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = c(0.1, 0.2), pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector with length equal to the number of columns in X"
  )
})

test_that("multivariate_analysis_pipeline errors when maf is not numeric", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = rep("0.1", ncol(d$X)), pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector"
  )
})

test_that("multivariate_analysis_pipeline errors when maf values are below 0", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  bad_maf <- d$maf
  bad_maf[1] <- -0.1
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = bad_maf, pip_cutoff_to_skip = 0
    ),
    "maf values must be between 0 and 1"
  )
})

test_that("multivariate_analysis_pipeline errors when maf values exceed 1", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  bad_maf <- d$maf
  bad_maf[2] <- 1.5
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = bad_maf, pip_cutoff_to_skip = 0
    ),
    "maf values must be between 0 and 1"
  )
})

test_that("multivariate_analysis_pipeline accepts maf of all zeros past validation", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  zero_maf <- rep(0, ncol(d$X))
  # maf = 0 is between 0 and 1 so should pass maf validation;
  # the pipeline will proceed past validation into skip_conditions
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = zero_maf, pip_cutoff_to_skip = 0
    ),
    error = function(e) e
  )
  # If it errors, the error should NOT be about maf validation
  if (inherits(result, "error")) {
    expect_false(grepl("maf values must be between", result$message))
  }
})

test_that("multivariate_analysis_pipeline accepts maf of all ones past validation", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  one_maf <- rep(1, ncol(d$X))
  # maf = 1 is between 0 and 1 so passes maf validation
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = one_maf, pip_cutoff_to_skip = 0
    ),
    error = function(e) e
  )
  # If it errors, the error should NOT be about maf validation
  if (inherits(result, "error")) {
    expect_false(grepl("maf values must be between", result$message))
  }
})

# ===========================================================================
# mvsusieR dependency check
# ===========================================================================

test_that("multivariate_analysis_pipeline errors when mvsusieR is not installed", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "mvsusieR"
  )
})

test_that("multivariate_analysis_pipeline error message includes install URL", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "github.com/stephenslab/mvsusieR"
  )
})

# ===========================================================================
# pip_cutoff_to_skip validation (tested via pipeline)
#
# skip_conditions is a local function inside the pipeline, so it can only
# be reached through the pipeline. When mvsusieR is not installed, we can
# still test cases where skip_conditions runs before mvsusieR is needed
# (when pip_cutoff_to_skip = 0, skip_conditions keeps all columns and
# the pipeline proceeds to mvsusieR check).
# ===========================================================================

test_that("pipeline with pip_cutoff_to_skip=0 passes skip_conditions and reaches mvsusieR check", {
  skip_if(requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR is installed, skipping not-installed test")
  d <- make_mv_data()
  # pip_cutoff_to_skip = 0 means keep all columns (no susie call needed)
  # Should reach mvsusieR check
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "mvsusieR"
  )
})

test_that("pipeline with pip_cutoff_to_skip wrong length vector errors", {
  d <- make_mv_data() # Y has 3 columns
  # Provide a vector of length 2 (not 1 and not ncol(Y)=3)
  # This should fail in skip_conditions before mvsusieR is called
  # But skip_conditions is inside the pipeline after the mvsusieR check,
  # so if mvsusieR is not installed, it will error on mvsusieR first.
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = c(0.1, 0.2)
    ),
    "pip_cutoff_to_skip must be a single number or a vector of the same length"
  )
})

test_that("pipeline with pip_cutoff_to_skip vector matching ncol(Y) is accepted", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach skip_conditions")
  d <- make_mv_data() # Y has 3 columns
  # Provide a vector of length 3 with all zeros to bypass susie calls
  # This should pass skip_conditions validation
  # It may still fail later in the pipeline, but not on pip_cutoff_to_skip
  result <- tryCatch(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf,
      pip_cutoff_to_skip = rep(0, 3)
    ),
    error = function(e) e
  )
  # If it errors, the error should not be about pip_cutoff_to_skip
  if (inherits(result, "error")) {
    expect_false(grepl("pip_cutoff_to_skip", result$message))
  }
})

# ===========================================================================
# filter_X_Y_missing behavior (tested indirectly through pipeline)
#
# Since filter_X_Y_missing is a local function, we test its logic by
# providing data that triggers its filtering paths.
# ===========================================================================

test_that("pipeline returns empty list when Y has all-NA rows leaving no data", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot reach filter_X_Y_missing")
  n <- 10
  p <- 5
  r <- 2
  set.seed(55)
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  rownames(X) <- paste0("s", 1:n)
  colnames(X) <- paste0("chr1:", 1:p * 100, ":A:G")
  # All rows of Y are NA
  Y <- matrix(NA_real_, nrow = n, ncol = r)
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("c", 1:r)
  maf <- rep(0.3, p)
  result <- multivariate_analysis_pipeline(
    X = X, Y = Y, maf = maf, pip_cutoff_to_skip = 0
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

# ===========================================================================
# Edge case: single-column Y (should warn and return empty list)
# ===========================================================================

test_that("pipeline warns and returns empty list when Y has single column", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, cannot test single-column behavior")
  d <- make_mv_data(r = 1)
  # With pip_cutoff_to_skip = 0, skip_conditions keeps the single column,
  # then checks ncol(Y_filtered) <= 1 which triggers warning + NULL return
  expect_warning(
    result <- multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = d$maf, pip_cutoff_to_skip = 0
    ),
    "After filtering.*1 context left"
  )
  expect_true(is.list(result))
  expect_equal(length(result), 0)
})

test_that("multivariate_analysis_pipeline errors when maf is NULL", {
  skip_if(!requireNamespace("mvsusieR", quietly = TRUE),
          "mvsusieR not installed, input validation unreachable")
  d <- make_mv_data()
  expect_error(
    multivariate_analysis_pipeline(
      X = d$X, Y = d$Y, maf = NULL, pip_cutoff_to_skip = 0
    ),
    "maf must be a numeric vector"
  )
})
