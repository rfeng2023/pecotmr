context("encoloc_extended")

# ---- filter_and_order_coloc_results (internal) ----
test_that("filter_and_order_coloc_results orders by PP values", {
  coloc_results <- data.frame(
    snp = c("s1", "s2", "s3"),
    PP.H4.1 = c(0.1, 0.5, 0.4),
    PP.H4.2 = c(0.3, 0.2, 0.5)
  )
  result <- pecotmr:::filter_and_order_coloc_results(coloc_results)
  expect_length(result, 2)  # 2 credible sets
  # First set should be ordered by PP.H4.1 decreasingly
  expect_equal(result[[1]][1, 2], 0.5)
})

test_that("filter_and_order_coloc_results errors on single column", {
  coloc_results <- data.frame(snp = c("s1", "s2"))
  expect_error(pecotmr:::filter_and_order_coloc_results(coloc_results),
               "Insufficient number of columns")
})

# ---- calculate_cumsum (internal) ----
test_that("calculate_cumsum computes cumulative sums of second column", {
  df <- data.frame(snp = c("s1", "s2", "s3"), pp = c(0.5, 0.3, 0.2))
  result <- pecotmr:::calculate_cumsum(df)
  expect_equal(result, c(0.5, 0.8, 1.0))
})

# ---- coloc_post_processor ----
test_that("coloc_post_processor warns without LD_meta_file_path", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_warning(
    coloc_post_processor(coloc_res),
    "LD_meta_file_path not provided"
  )
})

test_that("coloc_post_processor errors with LD path but no region", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_error(
    coloc_post_processor(coloc_res, LD_meta_file_path = "/some/path"),
    "analysis_region is not provided"
  )
})

test_that("coloc_post_processor warns with region but no LD path", {
  coloc_res <- list(summary = data.frame(PP.H4.abf = 0.9))
  expect_warning(
    coloc_post_processor(coloc_res, analysis_region = "chr1:100-200"),
    "will not be used"
  )
})
