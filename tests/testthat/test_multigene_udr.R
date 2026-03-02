context("multigene_udr")

library(testthat)

# =========================================================================
# multigene_udr.R: multigene_udr (lines 14-88)
# =========================================================================

test_that("multigene_udr with NULL exclude_condition defaults to empty vector", {
  skip_if_not_installed("udr")
  skip_if_not_installed("mashr")

  # Mock all heavy dependencies
  local_mocked_bindings(
    load_multitrait_R_sumstat = function(...) {
      list(
        bhat = matrix(rnorm(10), nrow = 5, ncol = 2),
        sbhat = matrix(abs(rnorm(10, mean = 1)), nrow = 5, ncol = 2)
      )
    },
    mash_rand_null_sample = function(...) {
      list(
        random = list(
          bhat = matrix(rnorm(20), nrow = 10, ncol = 2),
          sbhat = matrix(abs(rnorm(20, mean = 1)), nrow = 10, ncol = 2)
        ),
        null = list(
          bhat = matrix(rnorm(10), nrow = 5, ncol = 2),
          sbhat = matrix(abs(rnorm(10, mean = 1)), nrow = 5, ncol = 2)
        )
      )
    },
    filter_invalid_summary_stat = function(mash_input, bhat, sbhat, btoz, ...) {
      # Add z-scores if btoz is TRUE
      b_key <- bhat
      s_key <- sbhat
      if (btoz) {
        z_key <- gsub("\\.b$", ".z", b_key)
        mash_input[[z_key]] <- mash_input[[b_key]] / mash_input[[s_key]]
      }
      return(mash_input)
    },
    mash_pipeline = function(mash_input, ...) {
      list(U = list(mock_prior = matrix(1, 2, 2)), w = c(mock_prior = 1))
    }
  )

  combined_susie_list <- list(
    extracted_regional_window_combined_susie_result = "mock",
    extracted_regional_window_combined_sumstats_result = "mock"
  )

  result <- multigene_udr(
    combined_susie_list = combined_susie_list,
    coverage = 0.95,
    independent_variant_list = list("var1", "var2"),
    n_random = 10,
    n_null = 5,
    seed = 42,
    exclude_condition = NULL
  )
  expect_true(!is.null(result))
})

test_that("multigene_udr returns NULL when strong.b has < 2 rows", {
  skip_if_not_installed("udr")
  skip_if_not_installed("mashr")

  # Mock dependencies returning small data
  local_mocked_bindings(
    load_multitrait_R_sumstat = function(...) {
      list(
        bhat = matrix(rnorm(2), nrow = 1, ncol = 2),
        sbhat = matrix(abs(rnorm(2, mean = 1)), nrow = 1, ncol = 2)
      )
    },
    mash_rand_null_sample = function(...) {
      list(
        random = list(
          bhat = matrix(rnorm(4), nrow = 2, ncol = 2),
          sbhat = matrix(abs(rnorm(4, mean = 1)), nrow = 2, ncol = 2)
        ),
        null = list(
          bhat = matrix(rnorm(4), nrow = 2, ncol = 2),
          sbhat = matrix(abs(rnorm(4, mean = 1)), nrow = 2, ncol = 2)
        )
      )
    },
    filter_invalid_summary_stat = function(mash_input, bhat, sbhat, btoz, ...) {
      b_key <- bhat
      s_key <- sbhat
      if (btoz) {
        z_key <- gsub("\\.b$", ".z", b_key)
        mash_input[[z_key]] <- mash_input[[b_key]] / mash_input[[s_key]]
      }
      return(mash_input)
    }
  )

  combined_susie_list <- list(
    extracted_regional_window_combined_susie_result = "mock",
    extracted_regional_window_combined_sumstats_result = "mock"
  )

  result <- multigene_udr(
    combined_susie_list = combined_susie_list,
    coverage = 0.95,
    independent_variant_list = list("var1"),
    n_random = 5,
    n_null = 3,
    seed = 42,
    exclude_condition = c("cond1")
  )
  expect_null(result)
})

test_that("multigene_udr returns NULL when strong.b has < 2 columns", {
  skip_if_not_installed("udr")
  skip_if_not_installed("mashr")

  local_mocked_bindings(
    load_multitrait_R_sumstat = function(...) {
      list(
        bhat = matrix(rnorm(5), nrow = 5, ncol = 1),
        sbhat = matrix(abs(rnorm(5, mean = 1)), nrow = 5, ncol = 1)
      )
    },
    mash_rand_null_sample = function(...) {
      list(
        random = list(
          bhat = matrix(rnorm(10), nrow = 10, ncol = 1),
          sbhat = matrix(abs(rnorm(10, mean = 1)), nrow = 10, ncol = 1)
        ),
        null = list(
          bhat = matrix(rnorm(5), nrow = 5, ncol = 1),
          sbhat = matrix(abs(rnorm(5, mean = 1)), nrow = 5, ncol = 1)
        )
      )
    },
    filter_invalid_summary_stat = function(mash_input, bhat, sbhat, btoz, ...) {
      b_key <- bhat
      s_key <- sbhat
      if (btoz) {
        z_key <- gsub("\\.b$", ".z", b_key)
        mash_input[[z_key]] <- mash_input[[b_key]] / mash_input[[s_key]]
      }
      return(mash_input)
    }
  )

  combined_susie_list <- list(
    extracted_regional_window_combined_susie_result = "mock",
    extracted_regional_window_combined_sumstats_result = "mock"
  )

  result <- multigene_udr(
    combined_susie_list = combined_susie_list,
    coverage = 0.95,
    independent_variant_list = list("var1"),
    n_random = 5,
    n_null = 3,
    seed = 42,
    exclude_condition = NULL
  )
  expect_null(result)
})

test_that("multigene_udr exclude_condition is passed through correctly", {
  skip_if_not_installed("udr")
  skip_if_not_installed("mashr")

  captured_exclude <- NULL
  local_mocked_bindings(
    load_multitrait_R_sumstat = function(..., exclude_condition = NULL) {
      captured_exclude <<- exclude_condition
      list(
        bhat = matrix(rnorm(6), nrow = 3, ncol = 2),
        sbhat = matrix(abs(rnorm(6, mean = 1)), nrow = 3, ncol = 2)
      )
    },
    mash_rand_null_sample = function(...) {
      list(
        random = list(
          bhat = matrix(rnorm(10), nrow = 5, ncol = 2),
          sbhat = matrix(abs(rnorm(10, mean = 1)), nrow = 5, ncol = 2)
        ),
        null = list(
          bhat = matrix(rnorm(6), nrow = 3, ncol = 2),
          sbhat = matrix(abs(rnorm(6, mean = 1)), nrow = 3, ncol = 2)
        )
      )
    },
    filter_invalid_summary_stat = function(mash_input, bhat, sbhat, btoz, ...) {
      b_key <- bhat
      s_key <- sbhat
      if (btoz) {
        z_key <- gsub("\\.b$", ".z", b_key)
        mash_input[[z_key]] <- mash_input[[b_key]] / mash_input[[s_key]]
      }
      return(mash_input)
    },
    mash_pipeline = function(...) {
      list(U = list(mock = matrix(1, 2, 2)), w = c(mock = 1))
    }
  )

  combined_susie_list <- list(
    extracted_regional_window_combined_susie_result = "mock",
    extracted_regional_window_combined_sumstats_result = "mock"
  )

  result <- multigene_udr(
    combined_susie_list = combined_susie_list,
    coverage = 0.95,
    independent_variant_list = list("var1"),
    n_random = 5,
    n_null = 3,
    seed = 42,
    exclude_condition = c("brain", "liver")
  )
  expect_equal(captured_exclude, c("brain", "liver"))
})
