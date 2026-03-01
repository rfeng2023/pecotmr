context("quantile_twas_weight")

# ---- corr_filter (internal) ----
test_that("corr_filter removes highly correlated columns", {
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  # Make columns 1 and 2 highly correlated
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.01)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.9)
  expect_true(ncol(result$X.new) < p)
  expect_true(length(result$filter.id) == ncol(result$X.new))
})

test_that("corr_filter keeps all columns when uncorrelated", {
  set.seed(42)
  n <- 100; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  colnames(X) <- paste0("v", 1:p)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.99)
  expect_equal(ncol(result$X.new), p)
  expect_equal(result$filter.id, 1:p)
})

test_that("corr_filter preserves colnames for single remaining column", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n)
  colnames(X) <- c("a", "b", "c")
  # Make all columns nearly identical
  X[, 2] <- X[, 1] + rnorm(n, sd = 0.001)
  X[, 3] <- X[, 1] + rnorm(n, sd = 0.001)
  result <- pecotmr:::corr_filter(X, cor_thres = 0.5)
  expect_true(ncol(result$X.new) >= 1)
  expect_true(!is.null(colnames(result$X.new)))
})

# ---- calculate_coef_heterogeneity (internal) ----
test_that("calculate_coef_heterogeneity computes log(sd/mean)", {
  df <- data.frame(
    variant_id = c("v1", "v2"),
    coef_qr_0.25 = c(1, 0),
    coef_qr_0.50 = c(2, 0),
    coef_qr_0.75 = c(3, 0),
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 2)
  expect_true("coef_heter" %in% colnames(result))
  # v1: mean=2, sd=1, log(1/2) = -0.693
  expect_equal(result$coef_heter[1], log(1 / 2), tolerance = 0.01)
  # v2: mean=0, should be NA
  expect_true(is.na(result$coef_heter[2]))
})

test_that("calculate_coef_heterogeneity handles NA values", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = 1,
    coef_qr_0.50 = NA,
    coef_qr_0.75 = 3,
    stringsAsFactors = FALSE
  )
  result <- pecotmr:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 1)
  # Should handle NA gracefully
  expect_true(is.finite(result$coef_heter[1]) || is.na(result$coef_heter[1]))
})

# ---- calculate_xi_correlation ----
test_that("calculate_xi_correlation computes xi for valid monotonic data", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.05, 0.95, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- tau * 2  # monotonic increasing
  }
  result <- calculate_xi_correlation(df)
  expect_true("xi" %in% colnames(result))
  expect_true(is.numeric(result$xi))
})

test_that("calculate_xi_correlation warns on missing columns", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1", some_col = 1)
  expect_warning(result <- calculate_xi_correlation(df))
  expect_true(is.na(result$xi))
})

test_that("calculate_xi_correlation returns NA for too few valid values", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- NA
  }
  df$coef_qr_0.1 <- 1  # Only one valid value
  result <- calculate_xi_correlation(df, min_valid = 10)
  expect_true(is.na(result$xi))
})

# ---- qr_screen ----
test_that("qr_screen runs and returns results with quantreg", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  colnames(X) <- paste0("chr1:", seq(100, 500, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
})
