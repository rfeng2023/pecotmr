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
