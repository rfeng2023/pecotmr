context("dentist_qc")
library(MASS)
library(corpcor)

generate_dentist_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5, start_pos = 1000000, end_pos = 4000000) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    sumstat <- data.frame(
        position = unlist(lapply(seq(start_pos,end_pos,length.out = n_snps), round)),
        z = z_scores
    )
    return(list(sumstat = sumstat, LD_mat = ld_matrix, nSample = sample_size))
}

test_that("dentist output has exactly N rows for N input variants", {
    data <- generate_dentist_data(n_snps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentist output has exactly N rows with correct_chen_et_al_bug = FALSE", {
    data <- generate_dentist_data(n_snps = 100)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE))
    expect_equal(nrow(res), 100)
})

test_that("dentist output has exactly N rows for 200 variants", {
    data <- generate_dentist_data(seed = 99, n_snps = 200, sample_size = 200, n_outliers = 10)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 200)
})

test_that("dentist output has exactly N rows for 500 variants", {
    data <- generate_dentist_data(seed = 77, n_snps = 500, sample_size = 500, n_outliers = 25)
    expect_warning(res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 500)
})

test_that("dentist stops when missing position", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("something", "z")
    expect_error(dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE))
})

test_that("dentist stops when missing zscore", {
    data <- generate_dentist_data()
    colnames(data$sumstat) <- c("position", "something")
    expect_error(dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample, correct_chen_et_al_bug = FALSE))
})

generate_dentist_single_window_data <- function(seed=42, n_snps = 100, sample_size = 100, n_outliers = 5) {
    set.seed(seed)
    cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
    for (i in 1:(n_snps - 1)) {
        for (j in (i + 1):n_snps) {
            cor_matrix[i, j] <- runif(1, 0.2, 0.8)
            cor_matrix[j, i] <- cor_matrix[i, j]
        }
    }
    diag(cor_matrix) <- 1
    ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
    z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
    outlier_indices <- sample(1:n_snps, n_outliers)
    z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
    return(list(z_scores = z_scores, LD_mat = ld_matrix, nSample = sample_size))
}

test_that("dentist_single_window returns exactly N rows for N input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expect_equal(nrow(res), 100)
})

test_that("dentist_single_window warns when < 2000 variants", {
    data <- generate_dentist_single_window_data()
    expect_warning(dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
})

test_that("dentist_single_window returns exactly N rows with dedup", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample, duprThreshold = 0.99))
    expect_equal(nrow(res), 100)
})

test_that("dentist_single_window stops with zscore/LD matrix dimension mismatch", {
    data <- generate_dentist_single_window_data()
    expect_warning(expect_error(dentist_single_window(generate_dentist_single_window_data()$z_scores, R = generate_dentist_single_window_data(n_snps = 80)$LD_mat, nSample = data$nSample)))
})

test_that("dentist_single_window output columns are correct", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expected_cols <- c("original_z", "imputed_z", "iter_to_correct", "rsq", "is_duplicate", "outlier_stat", "outlier")
    expect_true(all(expected_cols %in% colnames(res)))
})

test_that("dentist_single_window original_z matches input z-scores", {
    data <- generate_dentist_single_window_data()
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expect_equal(res$original_z, data$z_scores)
})

test_that("add_dups_back_dentist works", {
    zScore <- c(1.2, 2.3, 2.4, 1.4, 5.6)
    dentist_output <- data.frame(
        original_z = c(1.2, 2.3, 5.6),
        imputed_z = c(1.1, 2.1, 5.1),
        iter_to_correct = c(1, 1, 3),
        rsq = c(0.9, 0.8,  0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )

    res <- pecotmr:::add_dups_back_dentist(zScore, dentist_output, find_dup_output)
    # Non-duplicates: z_diff is copied directly from dentist output
    expect_equal(res$z_diff[1], 0.1)
    expect_equal(res$z_diff[2], 0.2)
    expect_equal(res$z_diff[5], 0.5)
    # Duplicates: z_diff = (zScore - imputed_z) / sqrt(1 - rsq)
    expect_equal(res$z_diff[3], (2.4 - 2.1) / sqrt(1 - 0.8), tolerance = 1e-10)
    expect_equal(res$z_diff[4], (1.4 - (-1.1)) / sqrt(1 - 0.9), tolerance = 1e-10)
    expect_equal(res$imputed_z, c(1.1, 2.1, 2.1, -1.1, 5.1))
    expect_equal(res$is_duplicate, c(FALSE, FALSE, TRUE, TRUE, FALSE))
})

test_that("add_dups_back_dentist returns exactly N rows", {
    zScore <- c(1.2, 2.3, 2.4, 1.4, 5.6)
    dentist_output <- data.frame(
        original_z = c(1.2, 2.3, 5.6),
        imputed_z = c(1.1, 2.1, 5.1),
        iter_to_correct = c(1, 1, 3),
        rsq = c(0.9, 0.8, 0.5),
        z_diff = c(0.1, 0.2, 0.5)
    )
    find_dup_output <- data.frame(
        dupBearer = c(-1, -1, 2, 1, -1),
        sign = c(1, -1, 1, -1, 1)
    )
    res <- pecotmr:::add_dups_back_dentist(zScore, dentist_output, find_dup_output)
    expect_equal(nrow(res), 5)
})

test_that("add_dups_back_dentist stops when nrow mismatch", {
    z_scores <- rep(0, 5)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5, 6)
    )
    expect_error(pecotmr:::add_dups_back_dentist(z_scores, dentist_output, find_dup_output))
})

test_that("add_dups_back_dentist stops when zScore and find_dup_output nrow mismatch", {
    z_scores <- rep(0, 6)
    dentist_output <- list(
        original_z = c(1, 2, 3, 4, 5),
        imputed_z = c(1, 2, 3, 4, 5),
        iter_to_correct = c(1, 2, 3, 4, 5),
        rsq = c(1, 2, 3, 4, 5),
        z_diff = c(1, 2, 3, 4, 5)
    )
    find_dup_output <- list(
        dupBearer = c(-1, -1, -1, -1, -1),
        sign = c(1, 2, 3, 4, 5)
    )
    expect_error(pecotmr:::add_dups_back_dentist(z_scores, dentist_output, find_dup_output))
})

test_that("segment_by_dist works", {
    res <- pecotmr:::segment_by_dist(seq(2000000,5000000,100000), max_dist = 2000000, min_dim = 10)
    expect_true(nrow(res) >= 1)
})

test_that("segment_by_dist fill regions cover all input positions", {
    # Verify that fill regions cover every SNP index exactly once
    pos <- seq(1000000, 5000000, length.out = 200)
    res <- pecotmr:::segment_by_dist(pos, max_dist = 2000000, min_dim = 10)
    # Collect all fill region indices
    covered <- integer(0)
    for (k in 1:nrow(res)) {
        covered <- c(covered, res$fillStartIdx[k]:(res$fillEndIdx[k] - 1L))
    }
    # Every position from 1 to length(pos) should be covered
    expect_equal(sort(unique(covered)), 1:length(pos))
})

test_that("merge_windows returns exactly N rows", {
    data <- generate_dentist_data(n_snps = 1000, sample_size = 1000, start_pos = 0, end_pos = 2000)
    window_divided_res <- pecotmr:::segment_by_dist(data$sumstat$position, max_dist = 1000, min_dim = 10)
    dentist_result_by_window <- list()
    suppressWarnings({
        for (k in 1:nrow(window_divided_res)) {
            idx_range <- window_divided_res$windowStartIdx[k]:(window_divided_res$windowEndIdx[k] - 1L)
            zScore_k <- data$sumstat$z[idx_range]
            LD_mat_k <- data$LD_mat[idx_range, idx_range]
            dentist_result_by_window[[k]] <- dentist_single_window(
                zScore_k, R = LD_mat_k, nSample = 100,
                pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                ncpus = 1, correct_chen_et_al_bug = TRUE
            )
        }
    })
    res <- pecotmr:::merge_windows(dentist_result_by_window, window_divided_res)
    expect_equal(nrow(res), 1000)
})

test_that("merge_windows stops with window and imputed mismatch", {
    expect_error(pecotmr:::merge_windows(rep(0, 5), data.frame(windowStartIdx = rep(0,2), windowEndIdx = rep(0, 2))))
})

test_that("dentist windowed output has exactly N rows for large input", {
    # Generate data large enough to trigger windowed mode (> min_dim)
    data <- generate_dentist_data(seed = 123, n_snps = 1000, sample_size = 1000,
                                   n_outliers = 50, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                       min_dim = 100, window_size = 2000000)
    })
    expect_equal(nrow(res), 1000)
})

test_that("dentist windowed with correct_chen_et_al_bug=FALSE has exactly N rows", {
    data <- generate_dentist_data(seed = 456, n_snps = 1000, sample_size = 1000,
                                   n_outliers = 50, start_pos = 0, end_pos = 5000000)
    suppressWarnings({
        res <- dentist(data$sumstat, R = data$LD_mat, nSample = data$nSample,
                       min_dim = 100, window_size = 2000000, correct_chen_et_al_bug = FALSE)
    })
    expect_equal(nrow(res), 1000)
})

test_that("dentist outlier_stat formula is correct: (z-imputed)^2/(1-rsq)", {
    data <- generate_dentist_single_window_data(seed = 55, n_snps = 100)
    expect_warning(res <- dentist_single_window(data$z_scores, R = data$LD_mat, nSample = data$nSample))
    expected_stat <- (res$original_z - res$imputed_z)^2 / pmax(1 - res$rsq, 1e-8)
    expect_equal(res$outlier_stat, expected_stat, tolerance = 1e-10)
})

test_that("dentist with X matrix input returns exactly N rows", {
    set.seed(42)
    n_snps <- 80
    n_samples <- 100
    X <- matrix(rbinom(n_snps * n_samples, 2, 0.3), nrow = n_samples, ncol = n_snps)
    z_scores <- rnorm(n_snps)
    sumstat <- data.frame(position = seq(1000000, by = 1000, length.out = n_snps), z = z_scores)
    expect_warning(res <- dentist(sumstat, X = X))
    expect_equal(nrow(res), n_snps)
})

test_that("dentist_single_window with X matrix input returns exactly N rows", {
    set.seed(42)
    n_snps <- 80
    n_samples <- 100
    X <- matrix(rbinom(n_snps * n_samples, 2, 0.3), nrow = n_samples, ncol = n_snps)
    z_scores <- rnorm(n_snps)
    expect_warning(res <- dentist_single_window(z_scores, X = X))
    expect_equal(nrow(res), n_snps)
})

