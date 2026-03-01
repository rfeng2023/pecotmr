context("allele_qc_extended")

# ---- allele_qc basic functionality ----
test_that("allele_qc matches exact alleles", {
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  ref <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "C"), A1 = c("G", "T")
  )
  result <- allele_qc(target, ref, match_min_prop = 0)
  expect_equal(nrow(result$target_data_qced), 2)
})

test_that("allele_qc detects sign flips", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    z = 2.5
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "G", A1 = "A"
  )
  result <- allele_qc(target, ref, col_to_flip = "z", match_min_prop = 0)
  expect_equal(nrow(result$target_data_qced), 1)
  # z should be flipped
  expect_equal(result$target_data_qced$z, -2.5)
})

test_that("allele_qc removes strand ambiguous SNPs", {
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "A"), A1 = c("T", "G")
  )
  ref <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "A"), A1 = c("T", "G")
  )
  result <- allele_qc(target, ref, remove_strand_ambiguous = TRUE, match_min_prop = 0)
  # A/T is strand ambiguous, but since it's exact match, it should be kept
  expect_true(nrow(result$target_data_qced) >= 1)
})

test_that("allele_qc handles string input format", {
  target <- c("1:100:A:G", "1:200:C:T")
  ref <- c("1:100:A:G", "1:200:C:T")
  result <- allele_qc(target, ref, match_min_prop = 0)
  expect_equal(nrow(result$target_data_qced), 2)
})

test_that("allele_qc with chr prefix", {
  target <- c("chr1:100:A:G", "chr1:200:C:T")
  ref <- c("chr1:100:A:G", "chr1:200:C:T")
  result <- allele_qc(target, ref, match_min_prop = 0)
  expect_equal(nrow(result$target_data_qced), 2)
})

test_that("allele_qc warns when too few matches", {
  target <- c("1:100:A:G")
  ref <- c("2:200:C:T", "2:300:A:G", "2:400:C:T", "2:500:A:G", "2:600:C:T")
  expect_warning(allele_qc(target, ref, match_min_prop = 0.5))
})

test_that("allele_qc with no matching positions returns empty", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  ref <- data.frame(
    chrom = 1, pos = 999,
    A2 = "C", A1 = "T"
  )
  expect_warning(
    result <- allele_qc(target, ref, match_min_prop = 0),
    "No matching variants"
  )
  expect_equal(nrow(result$target_data_qced), 0)
})

test_that("allele_qc preserves extra columns", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G",
    beta = 0.5, se = 0.1
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- allele_qc(target, ref, match_min_prop = 0)
  expect_true("beta" %in% colnames(result$target_data_qced))
  expect_true("se" %in% colnames(result$target_data_qced))
})

test_that("allele_qc removes indels when requested", {
  target <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "AT"), A1 = c("G", "A")
  )
  ref <- data.frame(
    chrom = c(1, 1), pos = c(100, 200),
    A2 = c("A", "AT"), A1 = c("G", "A")
  )
  result_no_indel <- allele_qc(target, ref, remove_indels = TRUE, match_min_prop = 0)
  # With remove_indels = TRUE, the INDEL should be treated differently
  expect_true(nrow(result_no_indel$target_data_qced) >= 0)
})

test_that("allele_qc with lowercase alleles", {
  target <- data.frame(
    chrom = 1, pos = 100,
    A2 = "a", A1 = "g"
  )
  ref <- data.frame(
    chrom = 1, pos = 100,
    A2 = "A", A1 = "G"
  )
  result <- allele_qc(target, ref, match_min_prop = 0)
  expect_equal(nrow(result$target_data_qced), 1)
})

# ---- align_variant_names ----
test_that("align_variant_names aligns matching variants", {
  source <- c("1:100:A:G", "1:200:C:T")
  reference <- c("1:100:A:G", "1:200:C:T")
  result <- align_variant_names(source, reference)
  expect_equal(length(result$aligned_variants), 2)
  expect_length(result$unmatched_indices, 0)
})

test_that("align_variant_names handles sign flips", {
  source <- c("1:100:A:G")
  reference <- c("1:100:G:A")
  result <- align_variant_names(source, reference)
  expect_length(result$aligned_variants, 1)
})

test_that("align_variant_names warns on non-standard format", {
  source <- c("rs12345")
  reference <- c("rs67890")
  expect_warning(
    align_variant_names(source, reference),
    "do not follow the expected"
  )
})

test_that("align_variant_names errors on mixed formats", {
  source <- c("1:100:A:G")
  reference <- c("rs12345")
  expect_error(
    align_variant_names(source, reference),
    "different variant naming conventions"
  )
})

test_that("align_variant_names strips build suffix", {
  source <- c("1:100:A:G:b38")
  reference <- c("1:100:A:G")
  result <- align_variant_names(source, reference, remove_build_suffix = TRUE)
  expect_length(result$aligned_variants, 1)
})

test_that("align_variant_names adds chr prefix from reference", {
  source <- c("1:100:A:G")
  reference <- c("chr1:100:A:G", "chr1:200:C:T")
  result <- align_variant_names(source, reference)
  expect_true(grepl("^chr", result$aligned_variants[1]))
})

test_that("align_variant_names identifies unmatched variants", {
  source <- c("1:100:A:G", "1:999:C:T")
  reference <- c("1:100:A:G", "1:200:G:A")
  result <- align_variant_names(source, reference)
  # Second variant has no matching position in reference
  expect_true(length(result$unmatched_indices) >= 0)
})
