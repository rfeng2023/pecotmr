context("mash_wrapper_extended")

# ---- filter_invalid_summary_stat ----
test_that("filter_invalid_summary_stat replaces NaN/Inf in bhat", {
  dat <- list(
    bhat = data.frame(a = c(1, NaN, 3), b = c(Inf, 2, -Inf)),
    sbhat = data.frame(a = c(0.1, 0.2, 0.3), b = c(0.1, NA, 0.3))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat")
  expect_true(all(!is.nan(result$bhat)))
  expect_true(all(!is.infinite(result$bhat)))
  # NaN/Inf in bhat replaced with 0
  expect_equal(unname(result$bhat[1, 2]), 0)  # Inf -> 0
})

test_that("filter_invalid_summary_stat replaces NaN/Inf in sbhat", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.1, NaN, Inf))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat")
  # NaN/Inf in sbhat replaced with 1000
  expect_equal(unname(result$sbhat[1, "a"]), 0.1)
  expect_equal(unname(result$sbhat[2, "a"]), 1000)
  expect_equal(unname(result$sbhat[3, "a"]), 1000)
})

test_that("filter_invalid_summary_stat computes z from beta/se with btoz", {
  dat <- list(
    bhat = data.frame(a = c(1, 2, 3)),
    sbhat = data.frame(a = c(0.5, 1, 0.5))
  )
  result <- filter_invalid_summary_stat(dat, bhat = "bhat", sbhat = "sbhat", btoz = TRUE)
  expect_true("z" %in% names(result))
  expect_equal(as.numeric(result$z[1, 1]), 2.0)  # 1 / 0.5
})

test_that("filter_invalid_summary_stat processes z directly", {
  dat <- list(
    strong = list(z = data.frame(a = c(5, NaN, 0.1), b = c(1, 2, Inf))),
    random = list(z = data.frame(a = c(0.5, 0.2), b = c(0.3, 0.4)))
  )
  result <- filter_invalid_summary_stat(dat, z = "z")
  expect_true(all(!is.nan(result$strong$z)))
  expect_true(all(!is.infinite(result$strong$z)))
})

# ---- merge_mash_data ----
test_that("merge_mash_data combines two datasets", {
  d1 <- list(random = data.frame(a = 1:3, b = 4:6))
  d2 <- list(random = data.frame(a = 7:8, b = 9:10))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 5)
})

test_that("merge_mash_data handles NULL input", {
  d1 <- NULL
  d2 <- list(random = data.frame(a = 1:3))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("merge_mash_data handles empty input", {
  d1 <- list()
  d2 <- list(random = data.frame(a = 1:3))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

test_that("merge_mash_data handles different columns", {
  d1 <- list(random = data.frame(a = 1:2, b = 3:4))
  d2 <- list(random = data.frame(a = 5:6, c = 7:8))
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 4)
  expect_true(all(c("a", "b", "c") %in% colnames(result$random)))
})

test_that("merge_mash_data preserves data when one side empty", {
  d1 <- list(random = data.frame(a = 1:3))
  d2 <- list(random = NULL)
  result <- merge_mash_data(d1, d2)
  expect_equal(nrow(result$random), 3)
})

# ---- merge_susie_cs ----
test_that("merge_susie_cs handles single condition with single CS", {
  susie_fit <- list(list(
    cond1 = list(
      top_loci = data.frame(
        variant_id = c("1:100:A:G", "1:200:C:T"),
        pip = c(0.9, 0.1),
        cs_coverage_0.95 = c(1, 1),
        stringsAsFactors = FALSE
      )
    )
  ))

  result <- pecotmr:::merge_susie_cs(susie_fit)
  expect_s3_class(result, "data.frame")
  expect_true("variant_id" %in% colnames(result))
  expect_true("max_pip" %in% colnames(result))
})

# ---- extract_flatten_sumstats_from_nested ----
test_that("extract_flatten_sumstats_from_nested finds z scores", {
  data <- list(
    variant_names = c("1:100:A:G", "1:200:C:T"),
    sumstats = list(
      betahat = c(0.5, -0.3),
      sebetahat = c(0.1, 0.15)
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true("z" %in% colnames(result))
  expect_equal(result$z[1], 5.0)  # 0.5 / 0.1
})

test_that("extract_flatten_sumstats_from_nested extracts beta", {
  data <- list(
    variant_names = c("1:100:A:G"),
    sumstats = list(betahat = c(0.5), sebetahat = c(0.1))
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "beta")
  expect_equal(result$beta[1], 0.5)
})

test_that("extract_flatten_sumstats_from_nested extracts se", {
  data <- list(
    variant_names = c("1:100:A:G"),
    sumstats = list(betahat = c(0.5), sebetahat = c(0.1))
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "se")
  expect_equal(result$se[1], 0.1)
})

test_that("extract_flatten_sumstats_from_nested errors on invalid extract_inf", {
  data <- list()
  expect_error(extract_flatten_sumstats_from_nested(data, extract_inf = "invalid"),
               "must be one of")
})

test_that("extract_flatten_sumstats_from_nested adds chr prefix", {
  data <- list(
    variant_names = c("1:100:A:G"),
    sumstats = list(z = c(2.5))
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_true(grepl("^chr", result$variants[1]))
})

test_that("extract_flatten_sumstats_from_nested searches nested lists", {
  data <- list(
    level1 = list(
      variant_names = c("1:100:A:G"),
      sumstats = list(z = c(2.0))
    )
  )
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_s3_class(result, "data.frame")
  expect_equal(result$z[1], 2.0)
})

test_that("extract_flatten_sumstats_from_nested returns NULL when not found", {
  data <- list(a = list(b = 42))
  result <- extract_flatten_sumstats_from_nested(data, extract_inf = "z")
  expect_null(result)
})
