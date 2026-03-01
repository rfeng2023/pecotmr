context("load_phenotype_data region naming")

test_that("load_phenotype_data assigns colnames from region_name_col without extract_region_name", {
  # Create a temp phenotype file with a region name column
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5",
               "ENSG002\t300\t400\t3.5\t4.5"), tmp)

  # region_name_col = 1 means use the first column (gene_id) for naming
  result <- load_phenotype_data(tmp, region = NULL, region_name_col = 1)
  expect_type(result, "list")
  expect_length(result, 1)
  # After transposition, columns should be named with gene IDs
  expect_true(all(c("ENSG001", "ENSG002") %in% colnames(result[[1]])))
  file.remove(tmp)
})

test_that("load_phenotype_data works without region_name_col", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2",
               "ENSG001\t100\t200\t1.5\t2.5"), tmp)

  result <- load_phenotype_data(tmp, region = NULL)
  expect_type(result, "list")
  expect_length(result, 1)
  file.remove(tmp)
})

test_that("load_phenotype_data errors on empty phenotype file", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tstart\tend\tS1\tS2"), tmp)

  expect_error(
    load_phenotype_data(tmp, region = NULL),
    "empty"
  )
  file.remove(tmp)
})
