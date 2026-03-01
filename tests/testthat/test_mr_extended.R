context("mr_extended")

# ---- calc_I2 (internal) ----
test_that("calc_I2 returns 0 for small Q", {
  result <- pecotmr:::calc_I2(Q = list(0.0001), Est = c(1, 2, 3))
  expect_equal(result, 0)
})

test_that("calc_I2 returns valid value for large Q", {
  result <- pecotmr:::calc_I2(Q = list(10), Est = c(1, 2, 3))
  expect_true(result >= 0 && result <= 1)
})

test_that("calc_I2 clamps negative values to 0", {
  # Q just barely above threshold but (Q - Est + 1) / Q < 0
  result <- pecotmr:::calc_I2(Q = list(0.5), Est = c(1, 2, 3))
  expect_equal(result, 0)
})

# ---- mr_analysis ----
test_that("mr_analysis handles all-NA input", {
  input <- data.frame(
    gene_name = "GENE1",
    variant_id = NA_character_,
    bhat_x = NA_real_,
    sbhat_x = NA_real_,
    cs = NA_real_,
    pip = NA_real_,
    bhat_y = NA_real_,
    sbhat_y = NA_real_,
    stringsAsFactors = FALSE
  )

  result <- mr_analysis(input)
  expect_s3_class(result, "data.frame")
  expect_equal(result$gene_name, "GENE1")
  expect_true(is.na(result$meta_eff))
})

test_that("mr_analysis returns valid output for good input", {
  input <- data.frame(
    gene_name = rep("GENE1", 3),
    variant_id = c("1:100:A:G", "1:200:C:T", "1:300:G:A"),
    bhat_x = c(0.5, 0.3, 0.4),
    sbhat_x = c(0.1, 0.1, 0.1),
    cs = c(1, 1, 1),
    pip = c(0.4, 0.3, 0.3),
    bhat_y = c(0.2, 0.15, 0.18),
    sbhat_y = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mr_analysis(input, cpip_cutoff = 0.5)
  expect_s3_class(result, "data.frame")
  expect_true("meta_eff" %in% colnames(result))
  expect_true("meta_pval" %in% colnames(result))
  expect_true("Q" %in% colnames(result))
  expect_true("I2" %in% colnames(result))
  expect_equal(nrow(result), 1)
})

test_that("mr_analysis handles cpip below cutoff", {
  input <- data.frame(
    gene_name = "GENE1",
    variant_id = "1:100:A:G",
    bhat_x = 0.5,
    sbhat_x = 0.1,
    cs = 1,
    pip = 0.1,  # Low pip
    bhat_y = 0.2,
    sbhat_y = 0.05,
    stringsAsFactors = FALSE
  )

  result <- mr_analysis(input, cpip_cutoff = 0.5)
  # Should return null output since cpip < cutoff
  expect_true(is.na(result$meta_eff))
})

test_that("mr_analysis with multiple credible sets", {
  input <- data.frame(
    gene_name = rep("GENE1", 4),
    variant_id = paste0("1:", seq(100, 400, 100), ":A:G"),
    bhat_x = c(0.5, 0.3, 0.4, 0.6),
    sbhat_x = c(0.1, 0.1, 0.1, 0.1),
    cs = c(1, 1, 2, 2),
    pip = c(0.4, 0.6, 0.5, 0.5),
    bhat_y = c(0.2, 0.15, 0.3, 0.25),
    sbhat_y = c(0.05, 0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- mr_analysis(input, cpip_cutoff = 0.5)
  expect_equal(nrow(result), 1)
  expect_equal(result$num_CS, 2)
})
