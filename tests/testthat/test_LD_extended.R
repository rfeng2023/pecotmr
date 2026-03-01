context("LD_extended")

# ---- order_dedup_regions ----
test_that("order_dedup_regions orders and deduplicates", {
  df <- data.frame(
    chrom = c("chr2", "chr1", "chr1", "chr2"),
    start = c(100, 200, 100, 100),
    end = c(200, 300, 200, 200)
  )
  result <- pecotmr:::order_dedup_regions(df)
  expect_equal(nrow(result), 3)  # one duplicate removed
  expect_true(all(diff(result$start[result$chrom == result$chrom[1]]) >= 0))
})

test_that("order_dedup_regions strips chr prefix", {
  df <- data.frame(chrom = c("chr1", "chr2"), start = c(100, 200), end = c(200, 300))
  result <- pecotmr:::order_dedup_regions(df)
  expect_true(all(result$chrom %in% c(1L, 2L)))
})

# ---- find_intersection_rows ----
test_that("find_intersection_rows finds correct rows", {
  gd <- data.frame(
    chrom = c(1, 1, 1),
    start = c(0, 100, 200),
    end = c(100, 200, 300)
  )
  result <- pecotmr:::find_intersection_rows(gd, 1, 50, 250)
  expect_equal(result$start_row$start, 0)
  expect_equal(result$end_row$end, 300)
})

test_that("find_intersection_rows errors on uncovered region", {
  gd <- data.frame(
    chrom = c(1),
    start = c(100),
    end = c(200)
  )
  expect_error(pecotmr:::find_intersection_rows(gd, 1, 500, 600),
               "not covered")
})

test_that("find_intersection_rows adjusts bounds to available data", {
  gd <- data.frame(
    chrom = c(1),
    start = c(100),
    end = c(200)
  )
  # Region starts before available data
  result <- pecotmr:::find_intersection_rows(gd, 1, 50, 200)
  expect_equal(result$start_row$start, 100)
})

# ---- validate_selected_region ----
test_that("validate_selected_region passes for valid region", {
  start_row <- data.frame(start = 0)
  end_row <- data.frame(end = 300)
  expect_silent(pecotmr:::validate_selected_region(start_row, end_row, 50, 250))
})

test_that("validate_selected_region errors for uncovered region", {
  start_row <- data.frame(start = 100)
  end_row <- data.frame(end = 200)
  expect_error(
    pecotmr:::validate_selected_region(start_row, end_row, 50, 250),
    "not fully covered"
  )
})

# ---- extract_file_paths ----
test_that("extract_file_paths extracts correct paths", {
  gd <- data.frame(
    chrom = c(1, 1, 1),
    start = c(0, 100, 200),
    end = c(100, 200, 300),
    path = c("f1.ld", "f2.ld", "f3.ld")
  )
  intersection <- list(
    start_row = data.frame(chrom = 1, start = 0),
    end_row = data.frame(start = 200)
  )
  result <- pecotmr:::extract_file_paths(gd, intersection, "path")
  expect_equal(length(result), 3)
})

test_that("extract_file_paths errors on missing column", {
  gd <- data.frame(chrom = 1, start = 0, end = 100)
  intersection <- list(
    start_row = data.frame(chrom = 1, start = 0),
    end_row = data.frame(start = 0)
  )
  expect_error(pecotmr:::extract_file_paths(gd, intersection, "nonexistent"),
               "not found")
})

# ---- create_combined_LD_matrix ----
test_that("create_combined_LD_matrix merges non-overlapping blocks", {
  m1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("v1", "v2"), c("v1", "v2")))
  m2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2, dimnames = list(c("v3", "v4"), c("v3", "v4")))

  variants <- list(
    data.frame(variants = c("v1", "v2")),
    data.frame(variants = c("v3", "v4"))
  )
  result <- pecotmr:::create_combined_LD_matrix(list(m1, m2), variants)

  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_equal(result["v1", "v2"], 0.5)
  expect_equal(result["v3", "v4"], 0.3)
  expect_equal(result["v1", "v3"], 0)  # Cross-block should be 0
})

test_that("create_combined_LD_matrix handles overlapping boundary variant", {
  m1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(c("v1", "v2"), c("v1", "v2")))
  m2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2, dimnames = list(c("v2", "v3"), c("v2", "v3")))

  variants <- list(
    data.frame(variants = c("v1", "v2")),
    data.frame(variants = c("v2", "v3"))
  )
  result <- pecotmr:::create_combined_LD_matrix(list(m1, m2), variants)

  # v2 is shared, so total should be 3 variants
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

# ---- validate_block_structure ----
test_that("validate_block_structure passes for proper block structure", {
  # Create a proper block-diagonal matrix
  mat <- matrix(0, 6, 6)
  mat[1:3, 1:3] <- 0.5
  mat[4:6, 4:6] <- 0.5
  diag(mat) <- 1

  variant_ids <- paste0("v", 1:6)
  rownames(mat) <- colnames(mat) <- variant_ids

  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "1"),
    size = c(3, 3),
    start_idx = c(1, 4),
    end_idx = c(3, 6)
  )

  expect_silent(pecotmr:::validate_block_structure(mat, block_meta, variant_ids))
})

test_that("validate_block_structure errors on non-block structure", {
  mat <- matrix(0.5, 4, 4)
  diag(mat) <- 1

  variant_ids <- paste0("v", 1:4)
  rownames(mat) <- colnames(mat) <- variant_ids

  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "1"),
    size = c(2, 2),
    start_idx = c(1, 3),
    end_idx = c(2, 4)
  )

  expect_error(pecotmr:::validate_block_structure(mat, block_meta, variant_ids),
               "does not have the expected block structure")
})

# ---- merge_blocks (internal) ----
test_that("merge_blocks merges small adjacent blocks", {
  block_meta <- data.frame(
    block_id = c(1, 2, 3),
    chrom = c("1", "1", "1"),
    size = c(50, 50, 100),
    start_idx = c(1, 51, 101),
    end_idx = c(50, 100, 200)
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 100, max_size = 10000)
  expect_true(nrow(result) < 3)
})

test_that("merge_blocks does not merge cross-chromosome", {
  block_meta <- data.frame(
    block_id = c(1, 2),
    chrom = c("1", "2"),
    size = c(10, 10),
    start_idx = c(1, 11),
    end_idx = c(10, 20)
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 50, max_size = 10000)
  expect_equal(nrow(result), 2)  # Cannot merge across chromosomes
})

test_that("merge_blocks returns single block unchanged", {
  block_meta <- data.frame(
    block_id = 1, chrom = "1", size = 10,
    start_idx = 1, end_idx = 10
  )
  result <- pecotmr:::merge_blocks(block_meta, min_size = 100, max_size = 10000)
  expect_equal(nrow(result), 1)
})

# ---- can_merge (internal) ----
test_that("can_merge checks chromosome and size", {
  b1 <- data.frame(chrom = "1", size = 100)
  b2 <- data.frame(chrom = "1", size = 200)
  expect_true(pecotmr:::can_merge(b1, b2, max_size = 500))
  expect_false(pecotmr:::can_merge(b1, b2, max_size = 200))

  b3 <- data.frame(chrom = "2", size = 100)
  expect_false(pecotmr:::can_merge(b1, b3, max_size = 500))
})

# ---- partition_LD_matrix (internal) ----
test_that("partition_LD_matrix partitions correctly", {
  mat <- matrix(0, 6, 6)
  mat[1:3, 1:3] <- 0.5
  mat[4:6, 4:6] <- 0.5
  diag(mat) <- 1
  variant_ids <- paste0("v", 1:6)
  rownames(mat) <- colnames(mat) <- variant_ids

  ld_data <- list(
    combined_LD_matrix = mat,
    combined_LD_variants = variant_ids,
    block_metadata = data.frame(
      block_id = c(1, 2),
      chrom = c("1", "1"),
      size = c(3, 3),
      start_idx = c(1, 4),
      end_idx = c(3, 6)
    )
  )

  result <- pecotmr:::partition_LD_matrix(ld_data, merge_small_blocks = FALSE)

  expect_type(result, "list")
  expect_true("ld_matrices" %in% names(result))
  expect_true("variant_indices" %in% names(result))
  expect_length(result$ld_matrices, 2)
  expect_equal(nrow(result$ld_matrices[[1]]), 3)
  expect_equal(nrow(result$ld_matrices[[2]]), 3)
})

test_that("partition_LD_matrix errors on empty matrix", {
  ld_data <- list(
    combined_LD_matrix = matrix(nrow = 0, ncol = 0),
    combined_LD_variants = character(0),
    block_metadata = data.frame()
  )
  expect_error(pecotmr:::partition_LD_matrix(ld_data), "Empty or NULL")
})
