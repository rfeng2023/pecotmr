context("sumstats_qc")

test_that("rss_basic_qc requires correct columns", {
  sumstats <- data.frame(beta = 1, se = 0.5)
  LD_data <- list(combined_LD_variants = data.frame())
  expect_error(rss_basic_qc(sumstats, LD_data), "Missing columns")
})

test_that("rss_basic_qc processes matching variants correctly", {
  # Create matching sumstats and LD data
  variant_ids <- c("1:100:A:G", "1:200:C:T", "1:300:G:A")

  sumstats <- data.frame(
    chrom = c(1, 1, 1),
    pos = c(100, 200, 300),
    A1 = c("G", "T", "A"),
    A2 = c("A", "C", "G"),
    beta = c(0.5, -0.3, 0.1),
    se = c(0.1, 0.15, 0.2),
    z = c(5.0, -2.0, 0.5),
    stringsAsFactors = FALSE
  )

  LD_mat <- diag(3)
  rownames(LD_mat) <- colnames(LD_mat) <- variant_ids

  ref_panel <- data.frame(
    chrom = c(1, 1, 1),
    pos = c(100, 200, 300),
    A2 = c("A", "C", "G"),
    A1 = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  LD_data <- list(
    combined_LD_variants = ref_panel,
    combined_LD_matrix = LD_mat
  )

  result <- rss_basic_qc(sumstats, LD_data)
  expect_type(result, "list")
  expect_true("sumstats" %in% names(result))
  expect_true("LD_mat" %in% names(result))
})

test_that("summary_stats_qc errors on invalid method", {
  sumstats <- data.frame(variant_id = "1:100:A:G", z = 2.0)
  LD_data <- list(combined_LD_matrix = matrix(1, 1, 1, dimnames = list("1:100:A:G", "1:100:A:G")))
  expect_error(summary_stats_qc(sumstats, LD_data, method = "invalid"),
               "Invalid quality control method")
})
