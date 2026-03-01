context("LD")
library(tidyverse)

generate_dummy_data <- function() {
  region <- data.frame(
    chrom = "chr1",
    start = c(1000),
    end = c(1190)
  )
  meta_df <- data.frame(
    chrom = "chr1",
    start = c(1000, 1200, 1400, 1600, 1800),
    end = c(1200, 1400, 1600, 1800, 2000),
    path = c(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
      "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
      "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
      "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
      "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
    ))
  return(list(region = region, meta = meta_df))
}

# Generate a wider region that spans multiple blocks for partition testing
generate_multi_block_data <- function() {
  region <- data.frame(
    chrom = "chr1",
    start = c(1000),
    end = c(1500)
  )
  meta_df <- data.frame(
    chrom = "chr1",
    start = c(1000, 1200, 1400, 1600, 1800),
    end = c(1200, 1400, 1600, 1800, 2000),
    path = c(
      "./test_data/LD_block_1.chr1_1000_1200.float16.txt.xz,./test_data/LD_block_1.chr1_1000_1200.float16.bim",
      "./test_data/LD_block_2.chr1_1200_1400.float16.txt.xz,./test_data/LD_block_2.chr1_1200_1400.float16.bim",
      "./test_data/LD_block_3.chr1_1400_1600.float16.txt.xz,./test_data/LD_block_3.chr1_1400_1600.float16.bim",
      "./test_data/LD_block_4.chr1_1600_1800.float16.txt.xz,./test_data/LD_block_4.chr1_1600_1800.float16.bim",
      "./test_data/LD_block_5.chr1_1800_2000.float16.txt.xz,./test_data/LD_block_5.chr1_1800_2000.float16.bim"
    ))
  return(list(region = region, meta = meta_df))
}

test_that("Check that we correctly retrieve the names from the matrix",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  variants <- unlist(
    c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G", "chr1:1120:A:G", "chr1:1160:A:G"))
  expect_equal(
    unlist(res$combined_LD_variants),
    variants)
  file.remove(LD_meta_file_path)
})

test_that("Check that the LD block has the appropriate rownames and colnames",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  variants <- unlist(
    c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G", "chr1:1120:A:G", "chr1:1160:A:G"))
  expect_identical(rownames(res$combined_LD_matrix), variants)
  expect_identical(colnames(res$combined_LD_matrix), variants)
  file.remove(LD_meta_file_path)
})

test_that("Check that the LD block contains the correct information",{
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  res <- load_LD_matrix(LD_meta_file_path, region)
  # Variant names
  variants <- unlist(
    c("chr1:1000:A:G", "chr1:1040:A:G", "chr1:1080:A:G", "chr1:1120:A:G", "chr1:1160:A:G"))
  # Check LD Block 1
  ld_block_one <- res$combined_LD_matrix
  ld_block_one_original <- as.matrix(
    read_delim(
      "test_data/LD_block_1.chr1_1000_1200.float16.txt.xz",
      delim = " ", col_names = F))
  rownames(ld_block_one_original) <- colnames(ld_block_one_original) <- variants
  expect_equal(ld_block_one, ld_block_one_original)
  file.remove(LD_meta_file_path)
})

# New tests for partition_LD_matrix function

test_that("partition_LD_matrix correctly partitions a single block", {
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix first
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Now partition the matrix
  partitioned <- partition_LD_matrix(ld_data)
  
  # Expectations for single block case
  expect_equal(length(partitioned$ld_matrices), 1)
  expect_equal(nrow(partitioned$variant_indices), length(ld_data$combined_LD_variants))
  expect_equal(unique(partitioned$variant_indices$block_id), 1)
  expect_identical(rownames(partitioned$ld_matrices[[1]]), ld_data$combined_LD_variants)
  expect_identical(colnames(partitioned$ld_matrices[[1]]), ld_data$combined_LD_variants)
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix correctly partitions multiple blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Now partition the matrix without merging blocks
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)
  
  # Check if we have the correct number of blocks
  # Should have block 1 (1000-1200), block 2 (1200-1400), and block 3 (1400-1600)
  expected_block_count <- 3
  expect_equal(length(partitioned$ld_matrices), expected_block_count)
  
  # Check if all variants are assigned to blocks
  expect_equal(nrow(partitioned$variant_indices), length(ld_data$combined_LD_variants))
  
  # Check if block IDs are correct
  expect_setequal(unique(partitioned$variant_indices$block_id), 1:expected_block_count)
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix properly merges small blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Set min_merged_block_size high enough to force merging
  # Each test block likely has 5 variants (based on the existing test)
  min_block_size <- 10
  
  # Now partition the matrix with block merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = TRUE, 
                                    min_merged_block_size = min_block_size)
  
  # We expect fewer blocks after merging
  expect_lt(length(partitioned$ld_matrices), 3)
  
  # Check if all variants are still assigned to blocks
  expect_equal(nrow(partitioned$variant_indices), length(ld_data$combined_LD_variants))
  
  # Check if merged blocks are larger than min_block_size
  block_sizes <- sapply(partitioned$ld_matrices, nrow)
  expect_true(all(block_sizes >= min_block_size | block_sizes == length(ld_data$combined_LD_variants)))
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix respects max_merged_block_size", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Set max_merged_block_size to a small value to prevent merging all blocks
  # Each test block likely has 5 variants (based on the existing test)
  max_block_size <- 8
  
  # Now partition the matrix with restricted block size
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = TRUE,
                                    min_merged_block_size = 2,
                                    max_merged_block_size = max_block_size)
  
  # Check if no block exceeds max_block_size
  block_sizes <- sapply(partitioned$ld_matrices, nrow)
  expect_true(all(block_sizes <= max_block_size))
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix handles empty matrix gracefully", {
  # Create an empty LD data structure
  empty_ld_data <- list(
    combined_LD_matrix = matrix(0, 0, 0),
    combined_LD_variants = character(0),
    block_metadata = data.frame(
      block_id = integer(0),
      chrom = character(0),
      size = integer(0),
      start_idx = integer(0),
      end_idx = integer(0)
    )
  )
  
  # Expect an error for empty matrix
  expect_error(partition_LD_matrix(empty_ld_data), "Empty or NULL LD matrix provided")
})

test_that("partition_LD_matrix validates block structure properly", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix that spans multiple blocks
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Create an invalid block structure by modifying the block_metadata
  invalid_ld_data <- ld_data
  
  # Assuming we have at least 2 blocks:
  if(nrow(invalid_ld_data$block_metadata) >= 2) {
    # Create overlapping blocks with invalid start/end indices
    invalid_ld_data$block_metadata$start_idx[2] <- invalid_ld_data$block_metadata$start_idx[1]
    invalid_ld_data$block_metadata$end_idx[1] <- invalid_ld_data$block_metadata$end_idx[2]
    
    # Introduce non-zero elements between blocks to trigger validation error
    if(length(invalid_ld_data$combined_LD_variants) >= 2) {
      idx1 <- invalid_ld_data$block_metadata$start_idx[1]
      idx2 <- invalid_ld_data$block_metadata$start_idx[2] + 1
      if(idx1 <= length(invalid_ld_data$combined_LD_variants) && 
         idx2 <= length(invalid_ld_data$combined_LD_variants)) {
        var1 <- invalid_ld_data$combined_LD_variants[idx1]
        var2 <- invalid_ld_data$combined_LD_variants[idx2]
        invalid_ld_data$combined_LD_matrix[var1, var2] <- 0.5
      }
    }
    
    # Expect an error for invalid block structure
    expect_error(partition_LD_matrix(invalid_ld_data), "Matrix does not have the expected block structure")
  }
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix properly maps variants to blocks", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Partition without merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)
  
  # Check that each variant is mapped to the correct block
  for(i in seq_along(partitioned$ld_matrices)) {
    # Get variants in this block matrix
    block_variants <- rownames(partitioned$ld_matrices[[i]])
    
    # Find these variants in the variant_indices dataframe
    variant_block_ids <- partitioned$variant_indices$block_id[
      match(block_variants, partitioned$variant_indices$variant_id)]
    
    # All variants should be mapped to this block
    expect_true(all(variant_block_ids == i))
  }
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix handles row/column name mismatches", {
  data <- generate_dummy_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Create a version with mismatched rownames and colnames
  mismatched_ld_data <- ld_data
  rownames(mismatched_ld_data$combined_LD_matrix) <- NULL
  colnames(mismatched_ld_data$combined_LD_matrix) <- NULL
  
  # Should not error and should fix the names
  partitioned <- partition_LD_matrix(mismatched_ld_data)
  
  # Check if names are fixed
  expect_identical(rownames(partitioned$ld_matrices[[1]]), ld_data$combined_LD_variants)
  expect_identical(colnames(partitioned$ld_matrices[[1]]), ld_data$combined_LD_variants)
  
  file.remove(LD_meta_file_path)
})

test_that("partition_LD_matrix correctly extracts blocks based on metadata", {
  data <- generate_multi_block_data()
  region <- data$region
  LD_meta_file_path <- gsub("//", "/", tempfile(pattern = "ld_meta_file", tmpdir = tempdir(), fileext = ".RDS"))
  write_delim(data$meta, LD_meta_file_path, delim = "\t")
  
  # Load the LD matrix
  ld_data <- load_LD_matrix(LD_meta_file_path, region)
  
  # Partition without merging
  partitioned <- partition_LD_matrix(ld_data, merge_small_blocks = FALSE)
  
  # For each block, check if the extracted matrix matches the expected submatrix
  for(i in seq_along(partitioned$ld_matrices)) {
    block_info <- partitioned$block_metadata[i, ]
    start_idx <- block_info$start_idx
    end_idx <- block_info$end_idx
    
    # Skip if indices are invalid
    if(start_idx > length(ld_data$combined_LD_variants) || 
       end_idx > length(ld_data$combined_LD_variants) ||
       end_idx < start_idx) next
    
    # Get variants for this block
    block_variants <- ld_data$combined_LD_variants[start_idx:end_idx]
    
    # Extract expected submatrix
    expected_submatrix <- ld_data$combined_LD_matrix[block_variants, block_variants, drop = FALSE]
    
    # Compare with actual block matrix
    expect_equal(partitioned$ld_matrices[[i]], expected_submatrix)
  }
  
  file.remove(LD_meta_file_path)
})


# Tests for order_dedup_regions function
test_that("order_dedup_regions removes duplicate regions", {
  # Create regions with duplicates
  regions_with_dups <- data.frame(
    chrom = c("chr1", "chr1", "chr1"),
    start = c(100, 100, 200),  # Note: first two rows are duplicates
    end = c(150, 150, 250)
  )
  
  result <- order_dedup_regions(regions_with_dups)
  # Should have removed duplicate and return only two rows
  expect_equal(nrow(result), 2)
  expect_equal(result$start, c(100, 200))
})

# Tests for find_intersection_rows function
test_that("find_intersection_rows correctly identifies start and end rows", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )
  
  # Region entirely within the dataset
  result <- find_intersection_rows(genomic_data, 1, 220, 330)
  expect_equal(result$start_row$start, 200)
  expect_equal(result$end_row$end, 350)
})

test_that("find_intersection_rows adjusts region bounds if needed", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )
  
  # Region extends beyond the dataset
  result <- find_intersection_rows(genomic_data, 1, 50, 500)
  # Should adjust to the bounds of the dataset
  expect_equal(result$start_row$start, 100)
  expect_equal(result$end_row$end, 450)
})

test_that("find_intersection_rows errors for non-overlapping regions", {
  # Create a simple genomic dataset
  genomic_data <- data.frame(
    chrom = c(1, 1, 1, 1),
    start = c(100, 200, 300, 400),
    end = c(150, 250, 350, 450)
  )
  
  # Region entirely outside the dataset
  expect_error(
    find_intersection_rows(genomic_data, 2, 100, 200),
    "Region of interest is not covered by any rows in the data frame."
  )
})

# Additional edge case tests for partition_LD_matrix
test_that("partition_LD_matrix handles blocks with different chromosomes", {
  # Create test data with blocks on different chromosomes
  test_matrix <- matrix(0, 4, 4)
  diag(test_matrix) <- 1  # Set diagonal to 1
  rownames(test_matrix) <- colnames(test_matrix) <- c("1:100:A:G", "1:200:C:T", "2:100:G:A", "2:200:T:C")
  
  block_metadata <- data.frame(
    block_id = 1:2,
    chrom = c(1, 2),
    size = c(2, 2),
    start_idx = c(1, 3),
    end_idx = c(2, 4)
  )
  
  test_ld_data <- list(
    combined_LD_matrix = test_matrix,
    combined_LD_variants = c("1:100:A:G", "1:200:C:T", "2:100:G:A", "2:200:T:C"),
    block_metadata = block_metadata
  )
  
  # Partition the matrix
  partitioned <- partition_LD_matrix(test_ld_data)
  
  # Should not merge blocks from different chromosomes
  expect_equal(length(partitioned$ld_matrices), 2)
  
  # Each block should have the correct variants
  expect_equal(rownames(partitioned$ld_matrices[[1]]), c("1:100:A:G", "1:200:C:T"))
  expect_equal(rownames(partitioned$ld_matrices[[2]]), c("2:100:G:A", "2:200:T:C"))
})

test_that("partition_LD_matrix works with edge case block structures", {
  # Test case: One large block and several tiny blocks that need merging
  large_block_size <- 15
  small_block_size <- 2
  
  # Create a matrix with blocks of varying sizes
  n_variants <- large_block_size + small_block_size * 3
  test_matrix <- matrix(0, n_variants, n_variants)
  # Set diagonal to 1
  diag(test_matrix) <- 1
  
  # Generate variant names
  variant_names <- paste0("1:", 100:(100+n_variants-1), ":A:G")
  rownames(test_matrix) <- colnames(test_matrix) <- variant_names
  
  # Create block metadata
  block_metadata <- data.frame(
    block_id = 1:4,
    chrom = rep(1, 4),
    size = c(large_block_size, small_block_size, small_block_size, small_block_size),
    start_idx = c(1, large_block_size+1, large_block_size+small_block_size+1, large_block_size+small_block_size*2+1),
    end_idx = c(large_block_size, large_block_size+small_block_size, large_block_size+small_block_size*2, n_variants)
  )
  
  test_ld_data <- list(
    combined_LD_matrix = test_matrix,
    combined_LD_variants = variant_names,
    block_metadata = block_metadata
  )
  
  # Set minimum block size to force merging of small blocks
  min_merged_size <- small_block_size + 1
  
  # Partition with merging
  partitioned <- partition_LD_matrix(test_ld_data, merge_small_blocks = TRUE, 
                                    min_merged_block_size = min_merged_size)
  
  # Should merge the small blocks but leave the large block alone
  expect_lt(length(partitioned$ld_matrices), 4)
  expect_gt(length(partitioned$ld_matrices), 1)
  
  # First block should still be large_block_size
  expect_equal(nrow(partitioned$ld_matrices[[1]]), large_block_size)
})

# Tests for extract_LD_for_region function
test_that("extract_LD_for_region extracts correct region", {
  # Create mock LD matrix and variants
  ld_variants <- data.frame(
    chrom = c(1, 1, 1, 1),
    variants = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    GD = NA,
    pos = c(100, 200, 300, 400),
    A1 = c("A", "C", "G", "T"),
    A2 = c("G", "T", "A", "C")
  )
  
  ld_matrix <- matrix(0, 4, 4)
  diag(ld_matrix) <- 1
  rownames(ld_matrix) <- colnames(ld_matrix) <- ld_variants$variants
  
  # Define a region that should include the middle two variants
  region <- data.frame(
    chrom = 1,
    start = 180,
    end = 320
  )
  
  result <- extract_LD_for_region(ld_matrix, ld_variants, region, NULL)
  
  # Should have extracted only the relevant variants
  expect_equal(nrow(result$extracted_LD_variants), 2)
  expect_equal(result$extracted_LD_variants$variants, c("1:200:C:T", "1:300:G:A"))
  
  # Matrix should be 2x2 with the correct row/column names
  expect_equal(dim(result$extracted_LD_matrix), c(2, 2))
  expect_equal(rownames(result$extracted_LD_matrix), c("1:200:C:T", "1:300:G:A"))
})

test_that("extract_LD_for_region works with extract_coordinates", {
  # Create mock LD matrix and variants
  ld_variants <- data.frame(
    chrom = c(1, 1, 1, 1),
    variants = c("1:100:A:G", "1:200:C:T", "1:300:G:A", "1:400:T:C"),
    GD = NA,
    pos = c(100, 200, 300, 400),
    A1 = c("A", "C", "G", "T"),
    A2 = c("G", "T", "A", "C")
  )
  
  ld_matrix <- matrix(0, 4, 4)
  diag(ld_matrix) <- 1
  rownames(ld_matrix) <- colnames(ld_matrix) <- ld_variants$variants
  
  # Define a region that should include all variants
  region <- data.frame(
    chrom = 1,
    start = 50,
    end = 450
  )
  
  # Define specific coordinates to extract
  extract_coordinates <- data.frame(
    chrom = c(1, 1),
    pos = c(100, 300)
  )
  
  result <- extract_LD_for_region(ld_matrix, ld_variants, region, extract_coordinates)
  
  # Should have extracted only the specified coordinates
  expect_equal(nrow(result$extracted_LD_variants), 2)
  expect_equal(result$extracted_LD_variants$variants, c("1:100:A:G", "1:300:G:A"))
  
  # Matrix should be 2x2 with the correct row/column names
  expect_equal(dim(result$extracted_LD_matrix), c(2, 2))
  expect_equal(rownames(result$extracted_LD_matrix), c("1:100:A:G", "1:300:G:A"))
})

# Tests for create_combined_LD_matrix function
test_that("create_combined_LD_matrix correctly combines matrices", {
  # Create two simple LD matrices with some overlapping variants
  matrix1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  rownames(matrix1) <- colnames(matrix1) <- c("1:100:A:G", "1:200:C:T")
  
  matrix2 <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  rownames(matrix2) <- colnames(matrix2) <- c("1:200:C:T", "1:300:G:A")
  
  # Create variants lists
  variants1 <- data.frame(variants = c("1:100:A:G", "1:200:C:T"))
  variants2 <- data.frame(variants = c("1:200:C:T", "1:300:G:A"))
  
  # Combine matrices
  combined <- create_combined_LD_matrix(
    LD_matrices = list(matrix1, matrix2),
    variants = list(variants1, variants2)
  )
  
  # Should have created a 3x3 matrix with all unique variants
  expect_equal(dim(combined), c(3, 3))
  expect_equal(rownames(combined), c("1:100:A:G", "1:200:C:T", "1:300:G:A"))
  
  # Check that values from original matrices are preserved
  expect_equal(combined["1:100:A:G", "1:200:C:T"], 0.5)
  expect_equal(combined["1:200:C:T", "1:300:G:A"], 0.3)
  
  # Check diagonal values are 1 - without using diag() to avoid name issues
  expect_equal(combined[1,1], 1)
  expect_equal(combined[2,2], 1)
  expect_equal(combined[3,3], 1)
})

# Test for internal merge_blocks function via partition_LD_matrix
test_that("merge_blocks properly handles blocks at chromosome boundaries", {
  # Create test data with small blocks at chromosome boundaries
  test_matrix <- matrix(0, 6, 6)
  diag(test_matrix) <- 1
  variant_names <- c("1:900:A:G", "1:950:C:T", "2:100:G:A", "2:150:T:C", "3:100:A:G", "3:150:C:T")
  rownames(test_matrix) <- colnames(test_matrix) <- variant_names
  
  # Create block metadata with small blocks at chromosome boundaries
  block_metadata <- data.frame(
    block_id = 1:3,
    chrom = c(1, 2, 3),
    size = c(2, 2, 2),
    start_idx = c(1, 3, 5),
    end_idx = c(2, 4, 6)
  )
  
  test_ld_data <- list(
    combined_LD_matrix = test_matrix,
    combined_LD_variants = variant_names,
    block_metadata = block_metadata
  )
  
  # Set min block size to force merging attempts
  min_block_size <- 3
  
  # Partition with merging
  partitioned <- partition_LD_matrix(test_ld_data, merge_small_blocks = TRUE, 
                                     min_merged_block_size = min_block_size)
  
  # Should not merge blocks across chromosome boundaries
  expect_equal(length(partitioned$ld_matrices), 3)
  
  # Each block should match its chromosome
  for (i in 1:3) {
    block_variants <- rownames(partitioned$ld_matrices[[i]])
    chrom_from_variants <- unique(as.integer(sub(":.*", "", block_variants)))
    expect_equal(length(chrom_from_variants), 1)  # Should only have one chromosome per block
    expect_equal(chrom_from_variants, i)  # Should match the expected chromosome
  }
})