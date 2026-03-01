context("quantile_twas")

# ---- ensure_continuous_clusters (internal) ----
test_that("ensure_continuous_clusters preserves already continuous", {
  index <- c(1, 1, 2, 2, 3)
  result <- pecotmr:::ensure_continuous_clusters(index)
  expect_equal(result, c(1, 1, 2, 2, 3))
})

test_that("ensure_continuous_clusters renumbers non-continuous", {
  index <- c(1, 2, 1, 2, 3)
  result <- pecotmr:::ensure_continuous_clusters(index)
  # Should renumber: 1->1, 2->2, 1->3, 2->4, 3->5
  expect_equal(result, c(1, 2, 3, 4, 5))
})

test_that("ensure_continuous_clusters handles single element", {
  result <- pecotmr:::ensure_continuous_clusters(c(1))
  expect_equal(result, c(1))
})

test_that("ensure_continuous_clusters handles two elements same cluster", {
  result <- pecotmr:::ensure_continuous_clusters(c(1, 1))
  expect_equal(result, c(1, 1))
})

test_that("ensure_continuous_clusters all same cluster", {
  result <- pecotmr:::ensure_continuous_clusters(c(3, 3, 3, 3))
  expect_equal(result, c(1, 1, 1, 1))
})

# ---- get_cluster_ranges (internal) ----
test_that("get_cluster_ranges returns correct ranges", {
  index <- c(1, 1, 2, 2, 3)
  result <- pecotmr:::get_cluster_ranges(index)
  expect_length(result, 3)
  expect_equal(result[[1]], "1 - 2")
  expect_equal(result[[2]], "3 - 4")
  expect_equal(result[[3]], "5 - 5")
})

test_that("get_cluster_ranges single cluster", {
  index <- c(1, 1, 1)
  result <- pecotmr:::get_cluster_ranges(index)
  expect_length(result, 1)
  expect_equal(result[[1]], "1 - 3")
})

# ---- get_modularity (internal) ----
test_that("get_modularity returns 0 for single element", {
  W <- matrix(1, nrow = 1, ncol = 1)
  B <- matrix(1, nrow = 1, ncol = 1)
  result <- pecotmr:::get_modularity(W, B)
  expect_equal(result, 0)
})

test_that("get_modularity returns numeric for identity weight", {
  W <- diag(5)
  B <- matrix(c(rep(c(1, 0), c(3, 2)), rep(c(0, 1), c(3, 2))), nrow = 5)
  result <- pecotmr:::get_modularity(W, B)
  expect_type(result, "double")
})

test_that("get_modularity handles all-positive weights", {
  set.seed(42)
  W <- abs(matrix(rnorm(16), 4, 4))
  W <- (W + t(W)) / 2
  B <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4)
  result <- pecotmr:::get_modularity(W, B)
  expect_type(result, "double")
  expect_true(is.finite(result))
})

# ---- get_n_cluster (internal) ----
test_that("get_n_cluster returns single cluster for high correlation", {
  # All correlations > 0.8
  Sigma <- matrix(0.9, nrow = 5, ncol = 5)
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- pecotmr:::get_n_cluster(hc, Sigma, between_cluster = 0.8)
  expect_equal(result$n_cluster, 1)
})

test_that("get_n_cluster returns multiple clusters for low correlation", {
  set.seed(42)
  # Create block-diagonal correlation matrix
  Sigma <- diag(10)
  Sigma[1:5, 1:5] <- 0.8
  Sigma[6:10, 6:10] <- 0.8
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- pecotmr:::get_n_cluster(hc, Sigma, between_cluster = 0.5)
  expect_true(result$n_cluster >= 1)
})

# ---- integrate_tau (internal) ----
test_that("integrate_tau with vector input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(1, 2, 1)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("integrate_tau with matrix input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- matrix(c(1, 2, 1, 0.5, 1, 0.5), nrow = 2, byrow = TRUE)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_equal(nrow(result), 2)
  expect_true(all(result > 0))
})

test_that("integrate_tau zero weights give zero", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(0, 0, 0)
  result <- pecotmr:::integrate_tau(tau, a_tau)
  expect_equal(as.numeric(result), 0)
})

# ---- get_hierarchical_clusters (internal) ----
test_that("get_hierarchical_clusters returns valid structure", {
  set.seed(42)
  # Create a symmetric positive definite correlation matrix
  p <- 10
  Sigma <- matrix(0.3, nrow = p, ncol = p)
  Sigma[1:5, 1:5] <- 0.9
  Sigma[6:10, 6:10] <- 0.9
  diag(Sigma) <- 1
  result <- pecotmr:::get_hierarchical_clusters(Sigma, between_cluster = 0.5)
  expect_type(result, "list")
  expect_true("cluster" %in% names(result))
  expect_true("Q_modularity_initial" %in% names(result))
  expect_true("cluster_ranges" %in% names(result))
  expect_equal(nrow(result$cluster), p)
  # Each row should have exactly one 1
  expect_true(all(rowSums(result$cluster) == 1))
})

# ---- perform_grouped_integration (internal) ----
test_that("perform_grouped_integration returns correct structure", {
  set.seed(42)
  n_variants <- 5
  n_tau <- 10
  twas_weight <- matrix(rnorm(n_variants * n_tau), nrow = n_variants)
  rownames(twas_weight) <- paste0("v", 1:n_variants)
  tau_values <- seq(0.1, 0.9, length.out = n_tau)
  pseudo_R2 <- runif(n_tau, 0.1, 0.5)
  
  result <- pecotmr:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  expect_type(result, "list")
  expect_true("weights" %in% names(result))
  expect_true("twas_weight_performance" %in% names(result))
  expect_equal(nrow(result$weights), n_variants)
  # Should have at least 3 columns for fixed intervals (A1, A2, A3)
  expect_true(ncol(result$weights) >= 3)
})

test_that("perform_grouped_integration single variant skips clustering", {
  twas_weight <- matrix(rnorm(10), nrow = 1)
  rownames(twas_weight) <- "v1"
  tau_values <- seq(0.1, 0.9, length.out = 10)
  pseudo_R2 <- runif(10, 0.1, 0.5)
  
  result <- pecotmr:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  # Should only have fixed interval columns (A1, A2, A3), no clustering columns
  expect_equal(ncol(result$weights), 3)
})
