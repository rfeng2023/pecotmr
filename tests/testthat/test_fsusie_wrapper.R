context("fsusie_wrapper")

# ---- cal_purity ----
test_that("cal_purity with min method and single element CS", {
  set.seed(42)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_equal(result[[1]], 1)
})

test_that("cal_purity with min method and multi-element CS", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_length(result, 1)
  expect_true(result[[1]] >= 0 && result[[1]] <= 1)
})

test_that("cal_purity with non-min method returns three values", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2, 3))

  result <- pecotmr:::cal_purity(l_cs, X, method = "susie")
  expect_length(result[[1]], 3)  # min, mean, median
})

test_that("cal_purity with non-min method single element returns (1,1,1)", {
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  l_cs <- list(c(1))

  result <- pecotmr:::cal_purity(l_cs, X, method = "susie")
  expect_equal(result[[1]], c(1, 1, 1))
})

test_that("cal_purity with multiple credible sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  l_cs <- list(c(1, 2), c(5, 6, 7))

  result <- pecotmr:::cal_purity(l_cs, X, method = "min")
  expect_length(result, 2)
})

# ---- fsusie_get_cs ----
test_that("fsusie_get_cs creates susie-like sets", {
  set.seed(42)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)

  fSuSiE_obj <- list(
    cs = list(c(1, 2, 3), c(5, 6)),
    alpha = list(
      c(0.4, 0.3, 0.2, 0.05, 0.02, 0.01, 0.01, 0.005, 0.003, 0.002),
      c(0.01, 0.02, 0.02, 0.05, 0.45, 0.35, 0.05, 0.02, 0.02, 0.01)
    )
  )

  result <- fsusie_get_cs(fSuSiE_obj, X, requested_coverage = 0.95)

  expect_type(result, "list")
  expect_true("cs" %in% names(result))
  expect_true("purity" %in% names(result))
  expect_true("cs_index" %in% names(result))
  expect_true("coverage" %in% names(result))
  expect_true("requested_coverage" %in% names(result))
  expect_equal(result$requested_coverage, 0.95)
  expect_equal(length(result$cs), 2)
  expect_equal(names(result$cs), c("L1", "L2"))
})
