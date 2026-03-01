context("regularized_regression")

# ---- prs_cs ----
test_that("prs_cs errors on invalid LD input", {
  expect_error(prs_cs(bhat = rnorm(5), LD = "not_a_list", n = 100),
               "valid list of LD blocks")
})

test_that("prs_cs errors on non-positive sample size", {
  expect_error(prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
               "valid sample size")
})

test_that("prs_cs errors on mismatched maf length", {
  expect_error(
    prs_cs(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100, maf = rep(0.3, 3)),
    "same as 'maf'"
  )
})

test_that("prs_cs errors on mismatched bhat and LD dimensions", {
  expect_error(
    prs_cs(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the sum"
  )
})

test_that("prs_cs runs successfully with valid input", {
  set.seed(42)
  p <- 10
  n <- 100
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  for (i in 1:(p - 1)) {
    R[i, i + 1] <- 0.3
    R[i + 1, i] <- 0.3
  }
  result <- prs_cs(bhat = bhat, LD = list(blk1 = R), n = n,
                   maf = rep(0.3, p), n_iter = 50, n_burnin = 10, thin = 2)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
  expect_equal(length(result$beta_est), p)
})

# ---- sdpr ----
test_that("sdpr errors on mismatched bhat and LD dimensions", {
  expect_error(
    sdpr(bhat = rnorm(10), LD = list(blk1 = diag(5)), n = 100),
    "same as the length of bhat"
  )
})

test_that("sdpr errors on non-positive sample size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = -1),
    "positive integer"
  )
})

test_that("sdpr errors on invalid per_variant_sample_size", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         per_variant_sample_size = c(100, -1, 100, 100, 100)),
    "positive values"
  )
})

test_that("sdpr errors on invalid array values", {
  expect_error(
    sdpr(bhat = rnorm(5), LD = list(blk1 = diag(5)), n = 100,
         array = c(0, 1, 3, 1, 0)),
    "0, 1, or 2"
  )
})

test_that("sdpr runs successfully", {
  set.seed(42)
  p <- 10
  bhat <- rnorm(p, sd = 0.1)
  R <- diag(p)
  result <- sdpr(bhat = bhat, LD = list(blk1 = R), n = 100,
                 iter = 50, burn = 10, thin = 2, verbose = FALSE)
  expect_type(result, "list")
  expect_true("beta_est" %in% names(result))
})

# ---- susie_weights ----
test_that("susie_weights returns zeros when no alpha/mu/X_column_scale_factors", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_weights(susie_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

test_that("susie_weights errors on dimension mismatch", {
  mock_fit <- list(pip = c(0.1, 0.2))
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  expect_error(susie_weights(X = X, susie_fit = mock_fit), "Dimension mismatch")
})

# ---- susie_ash_weights ----
test_that("susie_ash_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.1, 0.2, 0.3))
  result <- susie_ash_weights(susie_ash_fit = mock_fit)
  expect_equal(result, rep(0, 3))
})

# ---- susie_inf_weights ----
test_that("susie_inf_weights returns zeros when no expected fields", {
  mock_fit <- list(pip = c(0.4, 0.5))
  result <- susie_inf_weights(susie_inf_fit = mock_fit)
  expect_equal(result, rep(0, 2))
})

# ---- glmnet_weights ----
test_that("glmnet_weights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnet_weights(X, y, alpha = 1)
  expect_equal(nrow(result), p)
  expect_equal(ncol(result), 1)
})

test_that("enet_weights computes elastic net weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- enet_weights(X, y)
  expect_equal(nrow(result), p)
})

test_that("lasso_weights computes LASSO weights", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- lasso_weights(X, y)
  expect_equal(nrow(result), p)
})

test_that("glmnet_weights handles zero-variance columns", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 50; p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 3] <- 5
  y <- X[, 1] * 0.5 + rnorm(n)
  result <- glmnet_weights(X, y, alpha = 1)
  expect_equal(result[3, 1], 0)
})

# ---- bayes_alphabet_weights ----
test_that("bayes_alphabet_weights errors on dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN"), "same number of rows")
})

test_that("bayes_alphabet_weights errors on covariate dimension mismatch", {
  skip_if_not_installed("qgg")
  X <- matrix(rnorm(100), nrow = 10)
  y <- rnorm(10)
  Z <- matrix(rnorm(15), nrow = 5)
  expect_error(bayes_alphabet_weights(X, y, method = "bayesN", Z = Z),
               "same number of rows")
})

# ---- gbayes_rss ----
test_that("gbayes_rss errors without LD matrix", {
  skip_if_not_installed("qgg")
  sumstats <- data.frame(beta = 0.1, se = 0.05, n = 1000)
  expect_error(gbayes_rss(sumstats = sumstats), "Must provide LD")
})

test_that("gbayes_rss errors on invalid method", {
  skip_if_not_installed("qgg")
  sumstats <- data.frame(beta = 0.1, se = 0.05, n = 1000)
  LD <- diag(1)
  expect_error(gbayes_rss(sumstats = sumstats, LD = LD, method = "invalidMethod"),
               "not valid")
})

test_that("gbayes_rss errors on non-dataframe input", {
  skip_if_not_installed("qgg")
  expect_error(gbayes_rss(sumstats = list(beta = 0.1), LD = diag(1)))
})

test_that("gbayes_rss errors on dimension mismatch", {
  skip_if_not_installed("qgg")
  sumstats <- data.frame(beta = c(0.1, 0.2), se = c(0.05, 0.1), n = c(1000, 1000))
  LD <- diag(3)
  expect_error(gbayes_rss(sumstats = sumstats, LD = LD), "must correspond")
})
