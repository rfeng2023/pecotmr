context("univariate_rss_diagnostics")

# ---- get_susie_result ----
test_that("get_susie_result returns NULL for empty input", {
  result <- get_susie_result(list())
  expect_null(result)
})

test_that("get_susie_result returns NULL when susie_result_trimmed missing", {
  result <- get_susie_result(list(some_data = 42))
  expect_null(result)
})

test_that("get_susie_result returns trimmed result when present", {
  mock_result <- list(pip = c(0.1, 0.5, 0.3), sets = list(cs = list()))
  con_data <- list(susie_result_trimmed = mock_result)
  result <- get_susie_result(con_data)
  expect_equal(result, mock_result)
})

# ---- extract_top_pip_info ----
test_that("extract_top_pip_info finds top PIP variant", {
  con_data <- list(
    susie_result_trimmed = list(pip = c(0.1, 0.7, 0.2)),
    variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
    sumstats = list(z = c(1.0, 3.5, -0.5))
  )
  result <- extract_top_pip_info(con_data)
  expect_equal(result$top_variant, "1:200:C:T")
  expect_equal(result$top_pip, 0.7)
  expect_equal(result$top_z, 3.5)
  expect_equal(result$top_variant_index, 2)
  expect_true(is.na(result$cs_name))
  expect_true(is.na(result$variants_per_cs))
})

test_that("extract_top_pip_info computes p_value from z", {
  con_data <- list(
    susie_result_trimmed = list(pip = c(0.9, 0.05, 0.05)),
    variant_names = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
    sumstats = list(z = c(5.0, 0.5, -0.3))
  )
  result <- extract_top_pip_info(con_data)
  expected_pval <- z_to_pvalue(5.0)
  expect_equal(result$p_value, expected_pval)
})

# ---- parse_cs_corr ----
test_that("parse_cs_corr handles NA correlations", {
  df <- data.frame(
    cs_name = "L1",
    top_pip = 0.9,
    cs_corr = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true("cs_corr_max" %in% colnames(result))
  expect_true("cs_corr_min" %in% colnames(result))
  expect_true(is.na(result$cs_corr_max))
  expect_true(is.na(result$cs_corr_min))
})

test_that("parse_cs_corr splits comma-separated correlations", {
  df <- data.frame(
    cs_name = c("L1", "L2"),
    top_pip = c(0.9, 0.3),
    cs_corr = c("1,0.3", "0.3,1"),
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true("cs_corr_1" %in% colnames(result))
  expect_true("cs_corr_2" %in% colnames(result))
  expect_equal(result$cs_corr_max[1], 0.3)
  expect_equal(result$cs_corr_min[1], 0.3)
})

test_that("parse_cs_corr handles empty string", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "",
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_true(is.na(result$cs_corr_max))
})

test_that("parse_cs_corr handles multiple correlations", {
  df <- data.frame(
    cs_name = "L1",
    cs_corr = "1,0.5,0.2",
    stringsAsFactors = FALSE
  )
  result <- parse_cs_corr(df)
  expect_equal(result$cs_corr_max, 0.5)
  expect_equal(result$cs_corr_min, 0.2)
  expect_equal(result$cs_corr_1, 1)
  expect_equal(result$cs_corr_2, 0.5)
  expect_equal(result$cs_corr_3, 0.2)
})
