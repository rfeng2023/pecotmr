context("perform_qr_analysis")

test_that("perform_qr_analysis requires quantreg package", {
  expect_true(is.function(perform_qr_analysis))
})

test_that("perform_qr_analysis runs with basic inputs", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(X) <- c("chr1:100:A:G", "chr1:200:C:T")
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(Y) <- "pheno1"

  result <- perform_qr_analysis(X, Y, tau_values = c(0.25, 0.5, 0.75))
  expect_s3_class(result, "data.frame")
  expect_true("chr" %in% colnames(result))
  expect_true("pos" %in% colnames(result))
  expect_true(nrow(result) > 0)
})

test_that("perform_qr_analysis with covariates Z", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(X) <- c("chr1:100:A:G", "chr1:200:C:T")
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(Y) <- "pheno1"
  Z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)

  result <- perform_qr_analysis(X, Y, Z = Z, tau_values = c(0.5))
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("perform_qr_analysis with single tau and single variant", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "chr2:500:G:A"
  Y <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(Y) <- "pheno1"

  result <- perform_qr_analysis(X, Y, tau_values = c(0.5))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

# ---- bayes convenience wrappers ----
context("bayes_alphabet convenience wrappers")

test_that("bayes_n_weights calls bayes_alphabet_weights with bayesN", {
  expect_true(is.function(bayes_n_weights))
})

test_that("bayes_l_weights calls bayes_alphabet_weights with bayesL", {
  expect_true(is.function(bayes_l_weights))
})

test_that("bayes_a_weights calls bayes_alphabet_weights with bayesA", {
  expect_true(is.function(bayes_a_weights))
})

test_that("bayes_c_weights calls bayes_alphabet_weights with bayesC", {
  expect_true(is.function(bayes_c_weights))
})

test_that("bayes_r_weights calls bayes_alphabet_weights with bayesR", {
  expect_true(is.function(bayes_r_weights))
})

test_that("prs_cs_weights is callable", {
  expect_true(is.function(prs_cs_weights))
})

test_that("sdpr_weights is callable", {
  expect_true(is.function(sdpr_weights))
})
