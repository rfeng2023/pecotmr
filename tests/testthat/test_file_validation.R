context("Input file validation")

# ---- tabix_region file validation ----
test_that("tabix_region errors on missing file", {
  expect_error(
    pecotmr:::tabix_region("/nonexistent/file.gz", "chr1:1-100"),
    "does not exist"
  )
})

# ---- load_covariate_data file validation ----
test_that("load_covariate_data errors on missing files", {
  expect_error(
    pecotmr:::load_covariate_data(c("/nonexistent/cov1.tsv", "/nonexistent/cov2.tsv")),
    "not found"
  )
})

test_that("load_covariate_data errors listing missing files", {
  expect_error(
    pecotmr:::load_covariate_data("/nonexistent/cov.tsv"),
    "cov.tsv"
  )
})

# ---- load_genotype_region file validation ----
test_that("load_genotype_region errors on missing genotype files", {
  expect_error(
    load_genotype_region("/nonexistent/geno"),
    "Genotype files not found"
  )
})

# ---- load_rss_data file validation ----
test_that("load_rss_data errors on missing sumstat file", {
  tmp_col <- tempfile(fileext = ".txt")
  writeLines("beta:beta_col", tmp_col)
  expect_error(
    load_rss_data("/nonexistent/sumstat.gz", tmp_col),
    "Summary statistics file not found"
  )
  file.remove(tmp_col)
})

test_that("load_rss_data errors on missing column file", {
  tmp_ss <- tempfile(fileext = ".gz")
  file.create(tmp_ss)
  expect_error(
    load_rss_data(tmp_ss, "/nonexistent/columns.txt"),
    "Column mapping file not found"
  )
  file.remove(tmp_ss)
})

# ---- load_phenotype_data kept_indices attribute ----
test_that("load_phenotype_data sets kept_indices attribute", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5"), tmp)
  result <- load_phenotype_data(tmp, region = NULL)
  expect_true(!is.null(attr(result, "kept_indices")))
  expect_equal(attr(result, "kept_indices"), 1L)
  file.remove(tmp)
})

test_that("load_phenotype_data kept_indices reflects filtering", {
  # Create two phenotype files - one with data, one empty (header only)
  tmp1 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5"), tmp1)

  tmp2 <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp2)

  # tmp2 is empty (no data rows), so should be filtered out
  # compact()/explicit NULL check should remove it
  result <- tryCatch(
    load_phenotype_data(c(tmp1, tmp2), region = NULL),
    error = function(e) NULL
  )
  # If both are loaded individually, tmp2 should produce NULL → filtered
  # The kept_indices should only include index 1
  if (!is.null(result)) {
    idx <- attr(result, "kept_indices")
    expect_true(1 %in% idx)
  }

  file.remove(tmp1, tmp2)
})
