context("slalom")

test_that("slalom basic output structure", {
  set.seed(42)
  n <- 50
  z <- rnorm(n)
  R <- diag(n)
  # Add some off-diagonal correlations
  for (i in 1:(n - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }

  result <- slalom(zScore = z, R = R)

  expect_type(result, "list")
  expect_named(result, c("data", "summary"))
  expect_s3_class(result$data, "data.frame")
  expect_true("original_z" %in% colnames(result$data))
  expect_true("prob" %in% colnames(result$data))
  expect_true("pvalue" %in% colnames(result$data))
  expect_true("outliers" %in% colnames(result$data))
  expect_true("nlog10p_dentist_s" %in% colnames(result$data))
  expect_equal(nrow(result$data), n)
})

test_that("slalom summary contains expected fields", {
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true("lead_pip_variant" %in% names(result$summary))
  expect_true("n_total" %in% names(result$summary))
  expect_true("n_r2" %in% names(result$summary))
  expect_true("n_dentist_s_outlier" %in% names(result$summary))
  expect_true("fraction" %in% names(result$summary))
  expect_true("max_pip" %in% names(result$summary))
  expect_true("cs_95" %in% names(result$summary))
  expect_true("cs_99" %in% names(result$summary))
  expect_equal(result$summary$n_total, n)
})

test_that("slalom probabilities sum to 1", {
  set.seed(42)
  n <- 30
  z <- rnorm(n)
  R <- diag(n)

  result <- slalom(zScore = z, R = R)
  expect_equal(sum(result$data$prob), 1, tolerance = 1e-10)
})

test_that("slalom with strong signal identifies lead variant", {
  set.seed(42)
  n <- 20
  z <- rnorm(n)
  z[5] <- -5  # Strong signal at position 5

  R <- diag(n)
  result <- slalom(zScore = z, R = R, lead_variant_choice = "pvalue")

  # With pvalue choice, lead should be variant 5 (lowest pnorm)
  expect_equal(result$summary$lead_pip_variant, 5)
})

test_that("slalom with ABF lead variant choice", {
  set.seed(42)
  n <- 10
  z <- rnorm(n, sd = 0.5)
  z[3] <- 4  # Strong signal
  R <- diag(n)

  result <- slalom(zScore = z, R = R, lead_variant_choice = "abf")

  # ABF should identify the variant with highest posterior probability
  expect_equal(result$summary$lead_pip_variant, which.max(result$data$prob))
})

test_that("slalom 95% CS is subset of 99% CS", {
  set.seed(42)
  n <- 15
  z <- rnorm(n)
  z[3] <- 3
  R <- diag(n)

  result <- slalom(zScore = z, R = R)

  expect_true(all(result$summary$cs_95 %in% result$summary$cs_99))
})

test_that("slalom errors on mismatched dimensions", {
  z <- rnorm(10)
  R <- diag(5)  # Wrong size
  expect_error(slalom(zScore = z, R = R))
})

test_that("slalom errors on non-square R", {
  z <- rnorm(10)
  R <- matrix(rnorm(50), nrow = 5, ncol = 10)
  expect_error(slalom(zScore = z, R = R))
})

test_that("slalom accepts X matrix instead of R", {
  set.seed(42)
  n_samples <- 100
  n_snps <- 10
  X <- matrix(sample(0:2, n_samples * n_snps, replace = TRUE), nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("snp", 1:n_snps)
  z <- rnorm(n_snps)

  result <- slalom(zScore = z, X = X)
  expect_type(result, "list")
  expect_equal(nrow(result$data), n_snps)
})

test_that("slalom with custom abf_prior_variance", {
  set.seed(42)
  n <- 10
  z <- rnorm(n)
  R <- diag(n)

  result1 <- slalom(zScore = z, R = R, abf_prior_variance = 0.04)
  result2 <- slalom(zScore = z, R = R, abf_prior_variance = 0.2)

  # Different priors should give different results
  expect_false(identical(result1$data$prob, result2$data$prob))
})

test_that("slalom outlier detection with correlated variants", {
  set.seed(42)
  n <- 10
  z <- c(3, 2.8, 0.1, -0.2, 0.3, -0.1, 0.2, 0.4, -0.3, 0.1)
  R <- diag(n)
  R[1, 2] <- R[2, 1] <- 0.9
  R[1, 3] <- R[3, 1] <- 0.05

  result <- slalom(zScore = z, R = R, r2_threshold = 0.6)
  # Should detect that correlated variants are consistent
  expect_type(result$data$outliers, "logical")
})
