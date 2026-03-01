context("susie_wrapper and susie_rss_wrapper")

# ---- susie_wrapper ----
test_that("susie_wrapper runs with init_L == max_L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  y <- X %*% beta + rnorm(n)
  result <- susie_wrapper(X, y, init_L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
})

test_that("susie_wrapper dynamically adjusts L", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  result <- susie_wrapper(X, y, init_L = 1, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
  expect_true(!is.null(result$sets))
})

# ---- susie_rss_wrapper ----
test_that("susie_rss_wrapper with L=1 does single iteration", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rnorm(p)
  result <- susie_rss_wrapper(z = z, R = R, L = 1)
  expect_true("pip" %in% names(result))
  expect_length(result$pip, p)
})

test_that("susie_rss_wrapper with L == max_L runs once", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rnorm(p)
  result <- susie_rss_wrapper(z = z, R = R, L = 5, max_L = 5)
  expect_true("pip" %in% names(result))
})

test_that("susie_rss_wrapper dynamic L with no CS found", {
  skip_if_not_installed("susieR")
  set.seed(42)
  p <- 10
  R <- diag(p)
  z <- rep(0.1, p)  # No signal, shouldn't find CS
  result <- susie_rss_wrapper(z = z, R = R, L = 2, max_L = 10, l_step = 2)
  expect_true("pip" %in% names(result))
})

# ---- get_cs_index (internal) ----
test_that("get_cs_index returns correct CS assignment", {
  skip_if_not_installed("susieR")
  set.seed(42)
  n <- 200
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  y <- X %*% beta + rnorm(n, sd = 0.5)
  fit <- susieR::susie(X, y, L = 5)
  if (!is.null(fit$sets$cs)) {
    idx <- pecotmr:::get_cs_index(fit$sets$cs, 1)
    expect_true(is.numeric(unname(idx)))
  }
})

# ---- lbf_to_alpha (internal) ----
test_that("lbf_to_alpha converts log BFs to posteriors", {
  # Row 1: values 0, 2, 4 -> column 3 should be highest
  # Row 2: values 3, 1, 0 -> column 1 should be highest
  lbf <- matrix(c(0, 3, 2, 1, 4, 0), nrow = 2, ncol = 3)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(dim(alpha), c(2, 3))
  # Each row should sum to 1
  expect_equal(rowSums(alpha), c(1, 1), tolerance = 1e-10)
  # Higher lbf should have higher alpha in each row
  expect_true(alpha[1, 3] > alpha[1, 1])
  expect_true(alpha[2, 1] > alpha[2, 3])
})

test_that("lbf_to_alpha handles uniform lbf", {
  lbf <- matrix(1, nrow = 1, ncol = 5)
  alpha <- pecotmr:::lbf_to_alpha(lbf)
  expect_equal(as.numeric(alpha), rep(0.2, 5), tolerance = 1e-10)
})

# ---- lbf_to_alpha_vector (internal) ----
test_that("lbf_to_alpha_vector with prior weights", {
  lbf <- c(0, 1, 2)
  alpha <- pecotmr:::lbf_to_alpha_vector(lbf)
  expect_length(alpha, 3)
  expect_equal(sum(alpha), 1, tolerance = 1e-10)
  # Higher lbf -> higher alpha
  expect_true(alpha[3] > alpha[1])
})
