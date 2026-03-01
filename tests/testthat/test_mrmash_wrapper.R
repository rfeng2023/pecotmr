context("mrmash_wrapper")

test_that("filter_mixture_components filters zero matrices", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat3 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.5, mat2 = 0.3, mat3 = 0.2)
  conditions_to_keep <- c("A", "B")

  result <- filter_mixture_components(conditions_to_keep, U, w)

  # mat2 should be removed (all zeros)
  expect_true(!"mat2" %in% names(result$U))
  # weights should be rescaled to maintain sum
  expect_equal(sum(result$w), sum(w), tolerance = 1e-10)
})

test_that("filter_mixture_components removes low weight components", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    mat2 = matrix(c(0.8, 0.3, 0.3, 0.9), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 0.999, mat2 = 0.00001)  # mat2 below default cutoff
  conditions_to_keep <- c("A", "B")

  result <- filter_mixture_components(conditions_to_keep, U, w, w_cutoff = 1e-04)
  expect_true(!"mat2" %in% names(result$U))
})

test_that("filter_mixture_components errors on missing condition", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  )
  w <- c(mat1 = 1.0)
  expect_error(filter_mixture_components(c("A", "C"), U, w), "not found in matrix")
})

test_that("filter_mixture_components subsets conditions", {
  U <- list(
    mat1 = matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3, 3,
                  dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  )
  w <- c(mat1 = 1.0)

  result <- filter_mixture_components(c("A", "B"), U, w)
  expect_equal(nrow(result$U[[1]]), 2)
  expect_equal(ncol(result$U[[1]]), 2)
})
