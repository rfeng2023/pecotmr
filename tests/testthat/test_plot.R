context("plot")

library(testthat)

# =========================================================================
# plot.R: venn function (lines 6-8)
# =========================================================================

test_that("venn function creates a plot object without error", {
  skip_if_not_installed("ggvenn")
  data <- list(
    SuSiE = c("gene1", "gene2", "gene3"),
    Lasso = c("gene2", "gene3", "gene4"),
    Enet = c("gene1", "gene4", "gene5"),
    MR.ASH = c("gene3", "gene5", "gene6")
  )
  # Use pdf(NULL) to suppress graphics device output
  pdf(NULL)
  result <- tryCatch(
    pecotmr:::venn(data),
    error = function(e) e
  )
  dev.off()
  # Check it returned something (not an error)
  expect_false(inherits(result, "error"))
})
