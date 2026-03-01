context("twas_predict and clean_context_names")

# ---- twas_predict ----
test_that("twas_predict multiplies X by weights", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  weights_list <- list(method1_weights = c(0.5, -0.5))
  result <- twas_predict(X, weights_list)
  expect_length(result, 1)
  expect_equal(names(result), "method1_predicted")
  expected <- X %*% c(0.5, -0.5)
  expect_equal(result[[1]], expected)
})

test_that("twas_predict handles multiple weight methods", {
  set.seed(42)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  weights_list <- list(
    lasso_weights = c(1, 0, -1),
    enet_weights = c(0.5, 0.3, 0.2),
    susie_weights = c(0, 0, 1)
  )
  result <- twas_predict(X, weights_list)
  expect_length(result, 3)
  expect_equal(names(result), c("lasso_predicted", "enet_predicted", "susie_predicted"))
  # Verify computation for one method
  expect_equal(result$lasso_predicted, X %*% c(1, 0, -1))
})

test_that("twas_predict with zero weights gives zero predictions", {
  X <- matrix(1:6, nrow = 2, ncol = 3)
  weights_list <- list(null_weights = c(0, 0, 0))
  result <- twas_predict(X, weights_list)
  expect_true(all(result$null_predicted == 0))
})

test_that("twas_predict with single variant", {
  X <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  weights_list <- list(single_weights = 2.0)
  result <- twas_predict(X, weights_list)
  expect_equal(as.numeric(result$single_predicted), c(2, 4, 6))
})

# ---- clean_context_names ----
test_that("clean_context_names removes gene suffix", {
  context <- c("brain_ENSG00000123", "liver_ENSG00000123")
  gene <- "ENSG00000123"
  result <- clean_context_names(context, gene)
  expect_equal(result, c("brain", "liver"))
})

test_that("clean_context_names handles multiple genes", {
  context <- c("brain_GENE1", "liver_GENE2", "heart_GENE1")
  gene <- c("GENE1", "GENE2")
  result <- clean_context_names(context, gene)
  expect_equal(result, c("brain", "liver", "heart"))
})

test_that("clean_context_names no match leaves context unchanged", {
  context <- c("brain_ABC", "liver_XYZ")
  gene <- "NONEXISTENT"
  result <- clean_context_names(context, gene)
  expect_equal(result, c("brain_ABC", "liver_XYZ"))
})

test_that("clean_context_names longer gene removed first", {
  context <- c("tissue_ENSG00000123456")
  gene <- c("ENSG00000123", "ENSG00000123456")
  result <- clean_context_names(context, gene)
  # Longer gene should be processed first
  expect_equal(result, "tissue")
})
