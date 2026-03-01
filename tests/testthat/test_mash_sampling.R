context("mash_rand_null_sample and merge_sumstats_matrices")

# ---- mash_rand_null_sample ----
test_that("mash_rand_null_sample with z scores returns random and null", {
  set.seed(42)
  dat <- list(
    z = data.frame(
      cond1 = c(5, 0.1, 0.2, 0.3, 0.5, 6, 0.1, 0.2, 0.4, 0.3),
      cond2 = c(0.2, 0.3, 0.1, 0.5, 0.4, 0.1, 0.3, 0.2, 0.1, 0.5)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 5, n_null = 3,
                                   exclude_condition = c(), seed = 123)
  expect_type(result, "list")
  expect_true("random" %in% names(result))
  expect_true("null" %in% names(result))
  expect_true("z" %in% names(result$random))
  expect_equal(nrow(result$random$z), 5)
})

test_that("mash_rand_null_sample with bhat/sbhat format", {
  set.seed(42)
  dat <- list(
    bhat = data.frame(
      cond1 = c(0.1, 0.05, 0.02, 0.01, 0.03),
      cond2 = c(0.02, 0.01, 0.03, 0.05, 0.04)
    ),
    sbhat = data.frame(
      cond1 = c(0.1, 0.1, 0.1, 0.1, 0.1),
      cond2 = c(0.1, 0.1, 0.1, 0.1, 0.1)
    )
  )
  result <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                   exclude_condition = c(), seed = 42)
  expect_type(result, "list")
  expect_true("random" %in% names(result))
  expect_true("bhat" %in% names(result$random))
  expect_true("sbhat" %in% names(result$random))
})

test_that("mash_rand_null_sample with seed is reproducible", {
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1)
    )
  )
  result1 <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                    exclude_condition = c(), seed = 42)
  result2 <- mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                                    exclude_condition = c(), seed = 42)
  expect_equal(result1$random$z, result2$random$z)
})

test_that("mash_rand_null_sample NULL input returns NULL", {
  result <- mash_rand_null_sample(NULL, n_random = 5, n_null = 3,
                                   exclude_condition = c())
  expect_null(result)
})

test_that("mash_rand_null_sample exclude_condition with column names errors due to numeric indexing", {
  # NOTE: The function checks names with %in% but uses -exclude_condition
  # which requires numeric indices. This is a known issue in the source.
  dat <- list(
    z = data.frame(
      cond1 = c(0.1, 0.2, 0.3, 0.4, 0.5),
      cond2 = c(0.5, 0.4, 0.3, 0.2, 0.1),
      cond3 = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
  )
  # Using string names triggers the bug (unary operator on string)
  expect_error(
    mash_rand_null_sample(dat, n_random = 3, n_null = 2,
                           exclude_condition = "cond3", seed = 42),
    "invalid argument to unary operator"
  )
})

# ---- merge_sumstats_matrices ----
test_that("merge_sumstats_matrices validates empty list", {
  expect_error(merge_sumstats_matrices(list(), "z"),
               "non-empty list")
})

test_that("merge_sumstats_matrices validates value_column type", {
  expect_error(merge_sumstats_matrices(list(data.frame(a = 1)), c("a", "b")),
               "single string")
})

test_that("merge_sumstats_matrices validates id_column type", {
  expect_error(merge_sumstats_matrices(list(data.frame(a = 1)), "a", id_column = 123),
               "single string")
})

test_that("merge_sumstats_matrices merges simple case", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  m2 <- data.frame(variants = c("chr1:100:A:G", "chr1:300:G:A"), z = c(2.0, 0.5))
  result <- merge_sumstats_matrices(list(m1, m2), value_column = "z")
  expect_s3_class(result, "data.frame")
  expect_true("variants" %in% colnames(result))
})

test_that("merge_sumstats_matrices single matrix returns valid result", {
  m1 <- data.frame(variants = c("chr1:100:A:G", "chr1:200:C:T"), z = c(1.5, -2.0))
  result <- merge_sumstats_matrices(list(m1), value_column = "z")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})
