#' Detect Outliers Using Dentist Algorithm
#'
#' DENTIST (Detecting Errors iN analyses of summary staTISTics) is a quality control
#' tool for GWAS summary data. It uses linkage disequilibrium (LD) information from a reference
#' panel to identify and correct problematic variants by comparing observed GWAS statistics to
#' predicted values. It can detect errors in genotyping/imputation, allelic errors, and
#' heterogeneity between GWAS and LD reference samples.
#'
#' @param sum_stat A data frame containing summary statistics, including 'pos' or 'position' and 'z' or 'zscore' columns.
#' @param LD_mat A matrix containing LD (linkage disequilibrium) information.
#' @param nSample The number of samples.
#' @param window_size The size of the window for dividing the genomic region. Default is 2000000.
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param seed The random seed for reproducibility. Default is 999.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#'
#' @return A data frame containing the imputed result and detected outliers.
#'
#' The returned data frame includes the following columns:
#'
#' \describe{
#'   \item{\code{original_z}}{The original z-score values from the input \code{sum_stat}.}
#'   \item{\code{imputed_z}}{The imputed z-score values computed by the Dentist algorithm.}
#'   \item{\code{rsq}}{The coefficient of determination (R-squared) between original and imputed z-scores.}
#'   \item{\code{iter_to_correct}}{The number of iterations required to correct the z-scores, if applicable.}
#'   \item{\code{index_within_window}}{The index of the observation within the window.}
#'   \item{\code{index_global}}{The global index of the observation.}
#'   \item{\code{outlier_stat}}{The computed statistical value based on the original and imputed z-scores and R-squared.}
#'   \item{\code{outlier}}{A logical indicator specifying whether the observation is identified as an outlier based on the statistical test.}
#' }
#'
#'
#' @examples
#' # Example usage of dentist
#' dentist(sum_stat, LD_mat, nSample)
#'
#' @details
#' # correct_chen_et_al_bug may affect the result in three parts compared with the original DENTIST method:
#' # 1. the way that window is divided
#'    when there is only one window and window size >= position range, it will run two windows
#'    where the first window is the whole range (e.g., 0-4122) and second is the subset of it, (e.g., 742-4122), which adds no information
#'    so if we set correct_chen_et_al_bug=TRUE and this is the scenario for window_size and positions, then it only runs on the whole window
#'    A more general description for this is:
#'        - if correct_chen_et_al_bug = FALSE, then there can never be a single window
#'        - if window_size is invalid (<=0 or >=range):
#'            - if correct_chen_et_al_bug==TRUE, run dentist single window, so just window (0-4122)
#'            - if correct_chen_et_al_bug==FALSE, run in the dentist original way, so window (0-4122) and (742-4122)
#'        - if window_size is valid:
#'            - if correct_chen_et_al_bug==TRUE, divide the window but don't add additional window if the first window covers all, so just (0-4122)
#'            - if correct_chen_et_al_bug==FALSE, divide the window in the dentist original way, so window (0-4122) and (742-4122)
#' 2. comparison between iteration index t and nIter (explained in the source code)
#' 3. !grouping_tmp (explained in the source code)
#'
#' @export
dentist <- function(sum_stat, LD_mat, nSample,
                    window_size = 2000000, pValueThreshold = 5.0369e-8, propSVD = 0.4, gcControl = FALSE,
                    nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99, ncpus = 1, seed = 999,
                    correct_chen_et_al_bug = TRUE, match_original = FALSE,
                    use_original_windowing = FALSE, min_dim = 2000) {
  # detect for column names and order by pos
  if (!any(tolower(c("pos", "position")) %in% tolower(colnames(sum_stat))) ||
    !any(tolower(c("z", "zscore")) %in% tolower(colnames(sum_stat)))) {
    stop("Input sum_stat is missing either 'pos'/'position' or 'z'/'zscore' column.")
  }
  # rename to common column name
  if (!tolower("pos") %in% tolower(colnames(sum_stat))) {
    colnames(sum_stat)[which(tolower(colnames(sum_stat)) %in% tolower(c("position")))] <- "pos"
  }

  if (!tolower("z") %in% tolower(colnames(sum_stat))) {
    colnames(sum_stat)[which(tolower(colnames(sum_stat)) %in% tolower(c("zscore")))] <- "z"
  }

  sum_stat <- sum_stat %>% arrange(pos)

  if (use_original_windowing) {
    # Use the original DENTIST C++ binary's segmentation algorithm
    window_divided_res <- segment_by_dist(sum_stat$pos, max_dist = window_size, min_dim = min_dim)
    dentist_result_by_window <- list()
    for (k in 1:nrow(window_divided_res)) {
      # windowEndIdx is 1-based exclusive (one past last element), so convert to
      # inclusive range by subtracting 1.
      idx_range <- window_divided_res$windowStartIdx[k]:(window_divided_res$windowEndIdx[k] - 1L)
      zScore_k <- sum_stat$z[idx_range]
      LD_mat_k <- LD_mat[idx_range, idx_range]
      dentist_result_by_window[[k]] <- dentist_single_window(
        zScore_k, LD_mat_k, nSample,
        pValueThreshold, propSVD, gcControl,
        nIter, gPvalueThreshold, duprThreshold,
        ncpus, seed, correct_chen_et_al_bug, match_original
      )
    }
    dentist_result <- merge_windows(dentist_result_by_window, window_divided_res)
  } else if (window_size <= 0 | ((window_size >= max(sum_stat$pos) - min(sum_stat$pos) | is.na(window_size)) & (correct_chen_et_al_bug == TRUE))) {
    dentist_result <- dentist_single_window(
      sum_stat$z, LD_mat, nSample,
      pValueThreshold, propSVD, gcControl,
      nIter, gPvalueThreshold, duprThreshold,
      ncpus, seed, correct_chen_et_al_bug, match_original
    )
  } else {
    # divide windows using the existing (non-original) method
    window_divided_res <- divide_into_windows(sum_stat$pos, window_size = window_size, correct_chen_et_al_bug = TRUE)
    # compute dentist result for each window
    dentist_result_by_window <- list()
    for (k in 1:nrow(window_divided_res)) {
      zScore_k <- sum_stat$z[window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]]
      LD_mat_k <- LD_mat[
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k],
        window_divided_res$windowStartIdx[k]:window_divided_res$windowEndIdx[k]
      ]
      dentist_result_by_window[[k]] <- dentist_single_window(
        zScore_k, LD_mat_k, nSample,
        pValueThreshold, propSVD, gcControl,
        nIter, gPvalueThreshold, duprThreshold,
        ncpus, seed, correct_chen_et_al_bug, match_original
      )
    }
    # merge single window result and generate a final dentist_result (similar to dentist_result above)
    dentist_result <- merge_windows(dentist_result_by_window, window_divided_res)
  }
  return(dentist_result)
}

#' Perform DENTIST on a single window
#'
#' This function performs imputation of summary statistics for a single genomic window
#' using the Dentist algorithm.
#'
#' @param zScore A numeric vector containing the z-score values for variants within the window.
#' @param LD_mat A square matrix containing linkage disequilibrium (LD) information for variants within the window.
#' @param nSample The total number of samples.
#' @param pValueThreshold The p-value threshold for significance. Default is 5e-8.
#' @param propSVD The proportion of singular value decomposition (SVD) to use. Default is 0.4.
#' @param gcControl Logical indicating whether genomic control should be applied. Default is FALSE.
#' @param nIter The number of iterations for the Dentist algorithm. Default is 10.
#' @param gPvalueThreshold The genomic p-value threshold for significance. Default is 0.05.
#' @param duprThreshold The absolute correlation r value threshold to be considered duplicate. Default is 0.99.
#' @param ncpus The number of CPU cores to use for parallel processing. Default is 1.
#' @param seed The random seed for reproducibility. Default is 999.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug. Default is TRUE.
#'
#' @return data frame includes columns representing the imputed summary statistics and outlier detected.
#'
#' @examples
#' # Example usage of dentist_impute_single_window
#' library(MASS)
#' library(corpcor)
#' set.seed(999)
#' # Set the number of SNPs, sample size, and number of outliers
#' n_snps <- 1000
#' sample_size <- 10000
#' n_outliers <- 5
#'
#' # Generate a correlation matrix with more off-diagonal correlation
#' cor_matrix <- matrix(0, nrow = n_snps, ncol = n_snps)
#' for (i in 1:(n_snps - 1)) {
#'   for (j in (i + 1):n_snps) {
#'     cor_matrix[i, j] <- runif(1, 0.2, 0.8) # Generate random correlations between 0.2 and 0.8
#'     cor_matrix[j, i] <- cor_matrix[i, j]
#'   }
#' }
#' diag(cor_matrix) <- 1
#'
#' # Convert the correlation matrix to a positive definite matrix
#' ld_matrix <- cov2cor(make.positive.definite(cor_matrix))
#'
#' # Simulate Z-scores based on the LD matrix
#' z_scores <- mvrnorm(n = 1, mu = rep(0, n_snps), Sigma = ld_matrix)
#'
#' # Introduce outliers
#' outlier_indices <- sample(1:n_snps, n_outliers)
#' z_scores[outlier_indices] <- rnorm(n_outliers, mean = 0, sd = 5)
#' dentist_single_window(zScore, LD_mat, nSample)
#'
#' @seealso
#' \code{\link{dentist}} for detecting outliers using the Dentist algorithm.
#'
#' @references
#' https://github.com/Yves-CHEN/DENTIST
#' @export
dentist_single_window <- function(zScore, LD_mat, nSample,
                                  pValueThreshold = 5e-8, propSVD = 0.4, gcControl = FALSE,
                                  nIter = 10, gPvalueThreshold = 0.05, duprThreshold = 0.99,
                                  ncpus = 1, seed = 999, correct_chen_et_al_bug = TRUE,
                                  match_original = FALSE) {
  calculate_stat <- function(impOp_zScores, impOp_imputed, impOp_rsq) {
    (impOp_zScores - impOp_imputed)^2 / (1 - impOp_rsq)
  }

  outlier_test <- function(stat, lambda, alpha = 5e-8) {
    minusLogPvalueChisq <- function(stat) {
      p <- pchisq(stat, df = 1, lower.tail = FALSE)
      return(-log10(p))
    }
    ifelse(minusLogPvalueChisq(stat / lambda) > -log10(alpha), TRUE, FALSE)
  }
  # Check that number of variants cannot be below 2000
  if (length(zScore) < 2000) {
    warning("The number of variants is below 2000. The algorithm may not work as expected, as suggested by the original DENTIST.")
  }
  # Check that LD_mat dimensions match the length of zScore
  if (!is.matrix(LD_mat) || nrow(LD_mat) != ncol(LD_mat) || nrow(LD_mat) != length(zScore)) {
    stop("LD_mat must be a square matrix with dimensions equal to the length of zScore.")
  }
  # Remove dups
  # The original DENTIST binary's dupThresh is an r-squared threshold (default 0.99).
  # It converts to a correlation threshold: rThreshold = round(sqrt(dupThresh)*1000)/1000
  # e.g. dupThresh=0.99 -> rThreshold=0.995
  # Our duprThreshold parameter follows the same convention (r-squared), so we convert here.
  org_Zscore <- zScore
  dedup_res <- NULL
  rThreshold <- round(sqrt(duprThreshold) * 1000) / 1000
  if (duprThreshold < 1.0) {
    dedup_res <- find_duplicate_variants(zScore, LD_mat, rThreshold)
    num_dup <- sum(dedup_res$dupBearer != -1)
    if (num_dup > 0) {
      message(paste(num_dup, "duplicated variants out of a total of", length(zScore), "were found at r threshold of", rThreshold))
    }
    zScore <- dedup_res$filteredZ
    LD_mat <- dedup_res$filteredLD
  }

  # Run the core C++ implementation, collecting warnings
  # The Rcpp code already caps rsq_eigen at 1.0 when it exceeds 1,
  # so the warning is informational. We collect and report them but do not error.
  rsq_warnings <- character(0)
  warning_handler <- function(w) {
    if (grepl("Adjusted rsq_eigen value exceeding 1", w$message)) {
      rsq_warnings <<- c(rsq_warnings, w$message)
      invokeRestart("muffleWarning")
    }
  }
  verbose_iter <- getOption("pecotmr.dentist.verbose", FALSE)
  res <- withCallingHandlers(
    dentist_iterative_impute(
      LD_mat, nSample, zScore,
      pValueThreshold, propSVD, gcControl, nIter,
      gPvalueThreshold, ncpus, seed, correct_chen_et_al_bug,
      match_original, verbose_iter
    ),
    warning = warning_handler
  )
  if (length(rsq_warnings) > 0) {
    warning(sprintf("%d rsq_eigen values exceeded 1 (capped at 1.0). Max reported: %s",
                    length(rsq_warnings), rsq_warnings[length(rsq_warnings)]))
  }
  res <- as.data.frame(res)
  # Recover dups
  if (duprThreshold < 1.0) {
    res <- add_dups_back_dentist(org_Zscore, res, dedup_res)
  }
  # detect outlier
  lambda_original <- 1
  res %>%
    mutate(
      outlier_stat = z_diff^2,
      outlier = outlier_test(outlier_stat, lambda_original)
    ) %>%
    select(-z_diff)
}

#' Add duplicates back to DENTIST output
#'
#' This function takes the output from the DENTIST algorithm and adds back the duplicated variants
#' based on the output from the `find_duplicate_variants` function.
#' @param zScore The original zScore
#' @param dentist_output A data frame containing the output from the DENTIST algorithm.
#' @param find_dup_output A list containing the output from the `find_duplicate_variants` function.
#'
#' @return A data frame with duplicated variants added back and an additional column indicating duplicates.
#'
#' @noRd
add_dups_back_dentist <- function(zScore, dentist_output, find_dup_output) {
  # Extract relevant columns from the DENTIST output
  original_z <- dentist_output$original_z
  imputed_z <- dentist_output$imputed_z
  iter_to_correct <- dentist_output$iter_to_correct
  rsq <- dentist_output$rsq
  z_diff <- dentist_output$z_diff

  # Extract output from find_duplicate_variants
  dupBearer <- find_dup_output$dupBearer
  sign <- find_dup_output$sign

  # Get the number of rows in dupBearer
  nrows_dup <- length(dupBearer)

  if (nrow(dentist_output) != sum(dupBearer == -1)) {
    stop("The number of rows in the input data does not match the occurrences of -1 in dupBearer.")
  }

  if (length(zScore) != nrows_dup) {
    stop("Input zScore and find_dup_output have inconsistent dimension")
  }

  # Initialize assignIdx vector
  count <- 1
  assignIdx <- rep(0, nrows_dup)

  for (i in seq_along(dupBearer)) {
    if (dupBearer[i] == -1) {
      assignIdx[i] <- count
      count <- count + 1
    } else {
      assignIdx[i] <- dupBearer[i]
    }
  }

  # Create a new data frame to store the updated values
  updated_data <- data.frame(
    original_z = numeric(nrows_dup),
    imputed_z = numeric(nrows_dup),
    iter_to_correct = numeric(nrows_dup),
    rsq = numeric(nrows_dup),
    z_diff = numeric(nrows_dup),
    is_duplicate = logical(nrows_dup)
  )

  for (i in seq_len(nrows_dup)) {
    updated_data$original_z[i] <- zScore[i]
    updated_data$iter_to_correct[i] <- iter_to_correct[assignIdx[i]]
    updated_data$rsq[i] <- rsq[assignIdx[i]]
    if (dupBearer[i] == -1) {
      # Non-duplicate: copy values directly from de-duplicated output
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]]
      updated_data$z_diff[i] <- z_diff[assignIdx[i]]
      updated_data$is_duplicate[i] <- FALSE
    } else {
      # Duplicate: sign-flip imputed_z and recompute z_diff from this SNP's own z-score.
      # The original binary computes output stat as (z - imputed)^2 / (1 - rsq) using each
      # SNP's own z-score (DENTIST.h line 706), not zScore_e^2 from the bearer. We must
      # recompute z_diff here so that z_diff^2 matches the binary's stat.
      updated_data$imputed_z[i] <- imputed_z[assignIdx[i]] * sign[i]
      denom <- sqrt(max(1 - updated_data$rsq[i], 1e-8))
      updated_data$z_diff[i] <- (zScore[i] - updated_data$imputed_z[i]) / denom
      updated_data$is_duplicate[i] <- TRUE
    }
  }

  return(updated_data)
}

#' Segment Genomic Region by Distance (Original DENTIST Algorithm)
#'
#' Implements the same windowing/segmentation algorithm as the original DENTIST C++ binary's
#' \code{segmentingByDist} function. Windows are created using quarter-distance SNP index
#' lookups, with gap detection for centromeres and large gaps.
#'
#' @param pos Integer vector of base pair positions (must be sorted).
#' @param max_dist Maximum distance (bp) between SNPs for windowing. Default is 2000000.
#' @param min_dim Minimum number of SNPs per window. Default is 2000.
#' @param verbose Logical; print segmentation info. Default is FALSE.
#'
#' @return A data frame with columns: windowIdx, windowStartIdx, windowEndIdx,
#'   fillStartIdx, fillEndIdx. Start indices are 1-based inclusive;
#'   end indices (windowEndIdx, fillEndIdx) are 1-based exclusive (one past last element),
#'   matching the C++ convention. Use \code{startIdx:(endIdx - 1)} for R inclusive ranges.
#'
#' @details
#' This is a faithful R translation of the C++ \code{segmentingByDist} function.
#' The algorithm:
#' \enumerate{
#'   \item Precomputes for each SNP: the index of the farthest SNP within \code{max_dist},
#'         and the index of the SNP at \code{max_dist/4} distance.
#'   \item Detects gaps > \code{max_dist/4} in the position vector (e.g., centromeres).
#'   \item Creates overlapping windows that slide by half the distance cutoff, with fill
#'         regions covering the inner three-quarters of each window.
#'   \item The first window's fill starts at the window start; the last window's fill
#'         ends at the window end.
#' }
#'
#' @seealso \code{\link{dentist_single_window}}, \code{\link{dentist}}
#'
#' @noRd
segment_by_dist <- function(pos, max_dist = 2000000, min_dim = 2000, verbose = FALSE) {
  n <- length(pos)
  if (n == 0) stop("No positions provided")

  cutoff <- max_dist
  minBlockSize <- 2000

  # Precompute nextIdx: for each SNP i, the farthest SNP index within cutoff distance.
  # C++ uses 0-based; we translate to 1-based. Key: loop boundaries must allow
  # j to reach n+1 (one past end) so that j-1 = n (last valid 1-based index).
  nextIdx <- integer(n)
  for (i in 1:n) {
    if (i == 1) {
      j <- 2
      while (j <= n && pos[j] - pos[1] < cutoff) j <- j + 1
      nextIdx[1] <- min(j, n)
    } else {
      j <- nextIdx[i - 1]
      while (j <= n && pos[j] - pos[i] < cutoff) j <- j + 1
      nextIdx[i] <- min(j, n)
    }
  }

  # Precompute quaterIdx: for each SNP i, the last SNP index within cutoff/4 distance.
  # C++ logic: starting from the previous quaterIdx value, advance j while
  # pos[j] < cutoff/4 + pos[i], then store j-1.
  quaterIdx <- integer(n)
  # First element: find largest index where pos < cutoff/4 + pos[1]
  j <- 1
  while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[1])) j <- j + 1
  quaterIdx[1] <- max(j - 1, 1L)
  # Rest: advance from previous value
  for (i in 2:n) {
    j <- quaterIdx[i - 1]
    while (j <= n && pos[j] < cutoff / 4 + as.numeric(pos[i])) j <- j + 1
    quaterIdx[i] <- max(j - 1, 1L)
  }
  # Clamp to valid range [1, n]
  quaterIdx <- pmin(quaterIdx, n)
  quaterIdx <- pmax(quaterIdx, 1L)

  # Helper to chain quaterIdx lookups (equivalent to quaterIdx[quaterIdx[x]] in C++)
  q1 <- function(x) quaterIdx[x]
  q2 <- function(x) quaterIdx[quaterIdx[x]]
  q3 <- function(x) quaterIdx[quaterIdx[quaterIdx[x]]]
  q4 <- function(x) quaterIdx[quaterIdx[quaterIdx[quaterIdx[x]]]]

  # Find gaps > cutoff/4
  diffs <- diff(pos)
  gapSizeThresh <- cutoff / 4
  allGaps <- c(1L)  # start of chromosome (1-based)
  for (i in seq_along(diffs)) {
    if (diffs[i] > gapSizeThresh) {
      allGaps <- c(allGaps, i + 1L)  # 1-based index of SNP after gap
    }
  }
  allGaps <- c(allGaps, n + 1L)  # end sentinel (one past last)

  if (verbose && length(allGaps) - 2 > 0) {
    message(sprintf("No. of gaps found: %d", length(allGaps) - 2))
    for (i in 2:(length(allGaps) - 1)) {
      message(sprintf("  Gap %d: %d - %d", i - 1, pos[allGaps[i] - 1], pos[allGaps[i]]))
    }
  }

  startList <- integer(0)
  endList <- integer(0)
  fillStartList <- integer(0)
  fillEndList <- integer(0)

  for (k in seq_len(length(allGaps) - 1)) {
    firstSegIdx <- length(startList) + 1  # 1-based
    rangeSize <- allGaps[k + 1] - allGaps[k]  # number of SNPs in this region
    if (rangeSize < minBlockSize / 2) next
    if (rangeSize - min_dim < 0) next

    startIdx <- allGaps[k]  # 1-based
    endIdx_region <- allGaps[k + 1]  # 1-based, one past last

    # endIdx for window = 4 quarter steps + 1 from startIdx
    endIdx <- min(q4(startIdx) + 1, endIdx_region)

    notLastInterval <- TRUE
    old_startIdx <- startIdx
    times <- 0

    repeat {
      times <- times + 1
      if (times > 400) stop("Windowing iteration limit exceeded")

      # Fill: quarter to three-quarters
      fillStartList <- c(fillStartList, q1(startIdx))
      # C++ uses exclusive end (i < fillEnd), so q3(startIdx) is already the
      # 1-based exclusive end. Keep it as-is for consistent convention.
      # The last segment's fillEnd will be reset below to endList[last].
      fillEndList <- c(fillEndList, q3(startIdx))

      if (endIdx_region <= endIdx) {
        notLastInterval <- FALSE
        # If last interval is small, startIdx goes back one step
        if (as.numeric(pos[min(endIdx - 1, n)]) - as.numeric(pos[q1(old_startIdx)]) < cutoff) {
          startIdx <- q1(old_startIdx)
        }
      }

      startList <- c(startList, startIdx)
      # Keep endIdx as 1-based exclusive (one-past-end), matching the C++ convention.
      # Cap at endIdx_region (also one-past-end) but do NOT cap at n — that would
      # convert the exclusive end to an inclusive end, creating mixed conventions.
      endList <- c(endList, min(endIdx, endIdx_region))

      # Update to next
      old_startIdx <- startIdx
      startIdx <- q2(startIdx)
      endIdx <- min(q4(startIdx) + 1, endIdx_region)

      if (!notLastInterval) break
    }

    # Reset first and last segment fill boundaries
    if (length(startList) >= firstSegIdx) {
      fillStartList[firstSegIdx] <- startList[firstSegIdx]
      fillEndList[length(fillEndList)] <- endList[length(endList)]
    }
  }

  if (length(startList) == 0) stop("No intervals created by segmentation")

  if (verbose) {
    message("Intervals:")
    for (i in seq_along(startList)) {
      message(sprintf("  %d: SNPs %d-%d (fill %d-%d)",
                      i, startList[i], endList[i], fillStartList[i], fillEndList[i]))
    }
  }

  # windowEndIdx and fillEndIdx are 1-based exclusive (one past last element).
  # Cap at n+1 (one past the last valid 1-based index).
  endList <- pmin(endList, n + 1L)
  fillEndList <- pmin(fillEndList, n + 1L)

  data.frame(
    windowIdx = seq_along(startList),
    windowStartIdx = startList,
    windowEndIdx = endList,
    fillStartIdx = fillStartList,
    fillEndIdx = fillEndList
  )
}

#' Divide Genomic Region into Windows
#'
#' This function divides a genomic region into windows based on the specified window size and other parameters.
#'
#' @param pos A numeric vector containing the positions of variants.
#' @param window_size The size of the window for dividing the genomic region.
#' @param correct_chen_et_al_bug Logical indicating whether to correct the Chen et al. bug.
#'
#' @return A data frame containing information about the divided windows, including start and end indices for both windows and fillers.
#'
#' @details
#' - If \code{correct_chen_et_al_bug = FALSE}, then there can never be a single window.
#' - If the \code{window_size} is invalid (<=0 or >=range):
#'   - If \code{correct_chen_et_al_bug==TRUE}, the function runs with a single window, covering the entire genomic region.
#'   - If \code{correct_chen_et_al_bug==FALSE}, the function runs in the original Dentist way, dividing the genomic region into multiple windows.
#' - If the \code{window_size} is valid:
#'   - If \code{correct_chen_et_al_bug==TRUE}, the function divides the window but doesn't add additional window if the first window covers all.
#'   - If \code{correct_chen_et_al_bug==FALSE}, the function divides the window in the original Dentist way.
#'
#' @seealso
#' \code{\link{dentist_single_window}} for detecting outlier for summary statistics for a single window using the Dentist algorithm.
#'
#' @noRd
divide_into_windows <- function(pos, window_size, correct_chen_et_al_bug) {
  windowStartIdx <- c()
  windowEndIdx <- c()
  fillStartIdx <- c()
  fillEndIdx <- c()
  input_pos_start <- min(pos)
  input_pos_end <- max(pos)
  pos_range <- input_pos_end - input_pos_start
  if (window_size <= 0 | (window_size >= pos_range & correct_chen_et_al_bug == TRUE)) {
    windowStartIdx <- 1
    windowEndIdx <- length(pos)
    fillStartIdx <- 1
    fillEndIdx <- length(pos)
  } else {
    for (i in 1:(2 * ceiling(pos_range / window_size))) {
      start_pointer <- ((i - 1) * window_size) * 0.5 + 1
      end_pointer <- start_pointer + window_size
      start_idx <- which.min(abs(pos - start_pointer - input_pos_start))
      windowStartIdx <- c(windowStartIdx, start_idx)
      end_idx <- which.min(abs(pos - end_pointer - input_pos_start))
      windowEndIdx <- c(windowEndIdx, end_idx)
    }
    # avoid the situation where the last window is too small
    cutoff <- which.min(abs(pos - input_pos_end + window_size))
    closest_idx <- which.min(abs(windowStartIdx - cutoff))
    windowStartIdx <- windowStartIdx[windowStartIdx <= cutoff | seq_along(windowStartIdx) == closest_idx]
    windowEndIdx <- windowEndIdx[1:length(windowStartIdx)]

    if (length(windowStartIdx) == 1 & (windowEndIdx[1] < 4122 | correct_chen_et_al_bug == FALSE)) {
      # to avoid only 1 window
      windowEndIdx <- c(windowEndIdx, length(pos))
      start_pointer <- (pos[windowEndIdx] - pos[windowStartIdx]) * 0.25 + pos[windowStartIdx] # position
      start_idx <- which.min(abs(pos - start_pointer))
      windowStartIdx <- c(windowStartIdx, start_idx)
    }
    # decide the fillStartIdx and fillEndIdx
    if (length(windowStartIdx) == 1) { # for single window (correct_chen_et_al_bug must be TRUE)
      fillStartIdx <- 1
      fillEndIdx <- length(pos)
    } else { # for multiple windows
      for (i in 1:length(windowStartIdx)) {
        if (i == 1 & i != length(windowStartIdx)) {
          # for the first window and not the last window, fill Start from 1 and fill End is 0.75 quantile
          fill_start_idx <- 1
          fillStartIdx <- c(fillStartIdx, fill_start_idx)
          fill_end_pointer <- 0.75 * window_size + input_pos_start
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        } else if (i != length(windowStartIdx)) {
          # for the windows in the middle, the start is the next variant of the previous fillendIdx, and the end is the 0.75 quantile
          fillStartIdx <- c(fillStartIdx, fillEndIdx[i - 1] + 1)
          fill_end_pointer <- pos[windowStartIdx[i]] + 0.75 * window_size
          fill_end_idx <- which.max(pos[pos <= fill_end_pointer])
          fillEndIdx <- c(fillEndIdx, fill_end_idx)
        }
      }
      # for the last window, the start is the next variant of the previous fillEndIdx, and the end is the last variant
      fillStartIdx <- c(fillStartIdx, fillEndIdx[length(fillEndIdx)] + 1)
      fillEndIdx <- c(fillEndIdx, length(pos))
    }
  }
  # combine window information into a data frame
  window_index <- 1:length(windowStartIdx)
  window_divided_res <- data.frame(
    windowIdx = window_index,
    windowStartIdx = windowStartIdx,
    windowEndIdx = windowEndIdx,
    fillStartIdx = fillStartIdx,
    fillEndIdx = fillEndIdx
  )
  # check if the divided windows are valid
  invalid_ends <- window_divided_res$fillStartIdx[1] != 1 | window_divided_res$fillEndIdx[nrow(window_divided_res)] != length(pos)
  if (nrow(window_divided_res) == 1 & invalid_ends) {
    stop("Invalid window divided!")
  } else if (nrow(window_divided_res) > 1) {
    invalid_middle <- 0
    for (i in 1:(nrow(window_divided_res) - 1)) {
      if (fillStartIdx[i + 1] - fillEndIdx[i] != 1) {
        invalid_middle <- invalid_middle + 1
      }
    }
    if (invalid_ends | invalid_middle != 0) {
      stop("Invalid window divided!")
    }
  }
  return(window_divided_res)
}

#' Merge dentist Results by Window
#'
#' This function merges DENTIST results by window into a single data frame.
#'
#' @param dentist_result_by_window A list containing imputed results for each window.
#' @param window_divided_res A data frame containing information about the divided windows.
#'
#' @return A data frame containing merged results.
#'
#' @details
#' The function checks if the number of imputed results matches the number of windows.
#' It then merges the results by window, adding an index within the window and a global index.
#' Finally, it extracts the results within the fillers and combines them into a single data frame.
#'
#' @noRd
merge_windows <- function(dentist_result_by_window, window_divided_res) {
  if (length(dentist_result_by_window) != nrow(window_divided_res)) {
    stop("Different number of windows and imputed results!")
  }
  merged_results <- c()
  for (k in 1:nrow(window_divided_res)) {
    imputed_k <- dentist_result_by_window[[k]]
    imputed_k$index_within_window <- seq(1:nrow(imputed_k))
    imputed_k <- imputed_k %>%
      mutate(index_global = index_within_window + window_divided_res$windowStartIdx[k] - 1)
    extracted_results <- imputed_k %>%
      filter(index_global >= window_divided_res$fillStartIdx[k] & index_global < window_divided_res$fillEndIdx[k])
    merged_results <- rbind(merged_results, extracted_results)
  }
  # Filter out un-imputed variants (imputed_z == 0 and rsq == 0) after merging,
  # matching the original DENTIST binary which skips these when writing output.
  merged_results <- merged_results %>%
    filter(!(imputed_z == 0 & rsq == 0))
  return(merged_results)
}

#' Read DENTIST-format Summary Statistics
#'
#' Reads a GWAS summary statistics file in the format expected by the DENTIST C++ binary
#' and returns it as an R data frame with z-scores computed.
#'
#' @param gwas_summary Path to summary statistics file (tab/space-separated,
#'   8 columns: SNP A1 A2 freq beta se p N). May be gzipped.
#' @return A data frame with columns: SNP, A1, A2, freq, beta, se, p, N, z.
#'
#' @export
read_dentist_sumstat <- function(gwas_summary) {
  if (!file.exists(gwas_summary)) {
    stop(paste0("Summary statistics file not found: ", gwas_summary))
  }
  # Handle gzipped files
  if (grepl("\\.gz$", gwas_summary)) {
    con <- gzfile(gwas_summary, "rt")
    on.exit(close(con), add = TRUE)
    ss <- read.table(con, header = TRUE, stringsAsFactors = FALSE)
  } else {
    ss <- read.table(gwas_summary, header = TRUE, stringsAsFactors = FALSE)
  }
  expected <- c("SNP", "A1", "A2", "freq", "beta", "se", "p", "N")
  # Case-insensitive matching
  matched <- sapply(expected, function(e) {
    idx <- which(tolower(colnames(ss)) == tolower(e))
    if (length(idx) == 0) return(NA_integer_)
    idx[1]
  })
  if (any(is.na(matched))) {
    missing <- expected[is.na(matched)]
    stop(paste0("Missing columns in summary file: ", paste(missing, collapse = ", "),
                ". Expected: ", paste(expected, collapse = ", ")))
  }
  ss <- ss[, matched, drop = FALSE]
  colnames(ss) <- expected
  ss$z <- ss$beta / ss$se
  return(ss)
}

#' Run R DENTIST Implementation on DENTIST-format Input Files
#'
#' Takes the same file-based inputs as the DENTIST C++ binary (summary stats file +
#' PLINK bfile) and runs the R \code{\link{dentist}} implementation. This allows
#' direct comparison between the R implementation and the C++ binary on identical data.
#'
#' This function reuses existing package utilities for file I/O, allele QC, and LD
#' computation: \code{\link{read_bim}} for PLINK bim files, \code{\link{allele_qc}}
#' for allele matching/flipping, \code{\link{load_genotype_region}} for genotype loading,
#' and \code{compute_LD} (from misc.R) for LD matrix computation with mean imputation
#' and Rfast::cora when available.
#'
#' @param gwas_summary Path to the GWAS summary statistics file (DENTIST format:
#'   tab-separated, 8 columns with header: SNP A1 A2 freq beta se p N). May be gzipped.
#' @param bfile PLINK binary file prefix (expects .bed/.bim/.fam files).
#' @param nSample Sample size. If NULL, taken from the N column of the summary file
#'   (using the median value). Default is NULL.
#' @param window_size Window size in base pairs. Default is 2000000.
#' @param pValueThreshold P-value threshold for outlier detection. Default is 5.0369e-8.
#' @param propSVD SVD truncation proportion. Default is 0.4.
#' @param gcControl Logical; apply genomic control. Default is FALSE.
#' @param nIter Number of QC iterations. Default is 10.
#' @param gPvalueThreshold GWAS p-value threshold for grouping. Default is 0.05.
#' @param duprThreshold LD r-squared threshold for duplicate detection. Default is 0.99.
#' @param ncpus Number of CPU threads. Default is 1.
#' @param seed Random seed. Default is 999.
#' @param correct_chen_et_al_bug Logical; correct known bugs in original DENTIST. Default is TRUE.
#' @param match_original Logical; use original DENTIST RNG and seed values. Default is FALSE.
#' @param use_original_windowing Logical; use original DENTIST C++ windowing algorithm. Default is FALSE.
#' @param min_dim Minimum number of SNPs per window (for original windowing). Default is 2000.
#' @param verbose Logical; print progress messages. Default is TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{result}{Data frame from \code{\link{dentist}} with outlier detection results.}
#'     \item{sum_stat}{Aligned summary statistics data frame (with SNP names and positions).}
#'     \item{LD_mat}{The LD correlation matrix used.}
#'   }
#'
#' @details
#' This function performs the full pipeline:
#' \enumerate{
#'   \item Reads the summary statistics file via \code{\link{read_dentist_sumstat}}.
#'   \item Reads the PLINK bim file via \code{read_bim} and matches SNPs by ID
#'         to obtain chromosome and position information.
#'   \item Aligns alleles using \code{\link{allele_qc}}, which handles strand flips,
#'         allele swaps, and z-score sign flipping.
#'   \item Loads genotype data via \code{\link{load_genotype_region}} and computes
#'         the LD correlation matrix via \code{compute_LD} (mean imputation + Rfast::cora).
#'   \item Calls \code{\link{dentist}} with the aligned data.
#' }
#'
#' The result includes the aligned summary statistics and LD matrix so they can be
#' reused for further analysis or debugging.
#'
#' @seealso \code{\link{dentist}}, \code{\link{read_dentist_sumstat}},
#'   \code{\link{allele_qc}}, \code{\link{parse_dentist_output}}
#'
#' @export
dentist_from_files <- function(gwas_summary,
                                bfile,
                                nSample = NULL,
                                window_size = 2000000,
                                pValueThreshold = 5.0369e-8,
                                propSVD = 0.4,
                                gcControl = FALSE,
                                nIter = 10,
                                gPvalueThreshold = 0.05,
                                duprThreshold = 0.99,
                                ncpus = 1,
                                seed = 999,
                                correct_chen_et_al_bug = TRUE,
                                match_original = FALSE,
                                use_original_windowing = FALSE,
                                min_dim = 2000,
                                verbose = TRUE) {
  # 1. Read summary stats
  if (verbose) message("Reading summary statistics...")
  ss <- read_dentist_sumstat(gwas_summary)
  if (verbose) message(sprintf("  %d variants read", nrow(ss)))

  # 2. Read bim using existing utility (read_bim expects .bed path)
  if (verbose) message("Reading reference panel bim...")
  bim <- read_bim(paste0(bfile, ".bed"))
  # bim columns: chrom, id, gpos, pos, a1, a0
  if (verbose) message(sprintf("  %d variants in reference panel", nrow(bim)))

  # 3. Match summary stats to bim by SNP ID to get chrom/pos
  common_ids <- intersect(ss$SNP, bim$id)
  if (length(common_ids) == 0) {
    stop("No common SNPs between summary statistics and reference panel.")
  }
  if (verbose) message(sprintf("  %d SNPs in common", length(common_ids)))

  ss_match <- ss[match(common_ids, ss$SNP), , drop = FALSE]
  bim_match <- bim[match(common_ids, bim$id), , drop = FALSE]

  # 4. Allele QC using existing allele_qc() utility
  # Build target df: sumstats alleles + chrom/pos from bim + extra columns to carry through
  target_df <- data.frame(
    chrom = as.integer(bim_match$chrom),
    pos = as.integer(bim_match$pos),
    A1 = ss_match$A1,
    A2 = ss_match$A2,
    z = ss_match$z,
    SNP = ss_match$SNP,
    N = ss_match$N,
    bim_id = bim_match$id,
    stringsAsFactors = FALSE
  )

  # Build ref df: bim alleles (a1 = effect allele, a0 = other allele)
  ref_df <- data.frame(
    chrom = as.integer(bim_match$chrom),
    pos = as.integer(bim_match$pos),
    A1 = bim_match$a1,
    A2 = bim_match$a0,
    stringsAsFactors = FALSE
  )

  if (verbose) message("Aligning alleles via allele_qc()...")
  qc_result <- allele_qc(
    target_data = target_df,
    ref_variants = ref_df,
    col_to_flip = "z",
    match_min_prop = 0,
    remove_dups = TRUE,
    remove_strand_ambiguous = TRUE
  )
  aligned <- qc_result$target_data_qced

  if (nrow(aligned) == 0) {
    stop("No variants remaining after allele QC.")
  }

  # Sort by position
  aligned <- aligned[order(aligned$pos), , drop = FALSE]
  rownames(aligned) <- NULL

  if (verbose) {
    message(sprintf("  %d variants after allele QC (from %d common)",
                    nrow(aligned), length(common_ids)))
  }

  # 6. Load genotypes and compute LD
  # When match_original = TRUE, load raw B-allele counts (as in the .bed file)
  # to match the binary's exact floating-point behavior. Otherwise use the
  # standard load_genotype_region() which returns A-allele counts (2 - B).
  if (verbose) message("Loading genotype data...")
  if (match_original) {
    if (!requireNamespace("snpStats", quietly = TRUE)) {
      stop("snpStats is required for match_original = TRUE")
    }
    geno <- snpStats::read.plink(bfile)
    X <- as(geno$genotypes, "numeric")  # B-allele counts, matching binary
  } else {
    X <- load_genotype_region(bfile, region = NULL)
  }

  # Determine sample size — use reference panel size (nrow of genotype matrix),
  # NOT the GWAS sample size from the summary stats N column.
  # The original DENTIST binary uses ref.N (number of samples in the .fam file)
  # for the K = min(idx.size(), nSample) * propSVD truncation.
  if (is.null(nSample)) {
    nSample <- as.integer(nrow(X))
    if (verbose) message(sprintf("Using reference panel sample size: N = %d", nSample))
  }

  # Subset to aligned SNPs using bim IDs
  snp_ids_for_geno <- aligned$bim_id
  missing_snps <- snp_ids_for_geno[!snp_ids_for_geno %in% colnames(X)]
  if (length(missing_snps) > 0) {
    warning(sprintf("%d aligned SNPs not found in genotype matrix; dropping them.",
                    length(missing_snps)))
    aligned <- aligned[snp_ids_for_geno %in% colnames(X), , drop = FALSE]
    snp_ids_for_geno <- aligned$bim_id
  }
  X_sub <- X[, snp_ids_for_geno, drop = FALSE]

  # Remove zero-variance columns (monomorphic SNPs) to avoid NaN correlations
  # Check variance BEFORE handling missing (using non-missing values)
  col_vars <- apply(X_sub, 2, function(x) {
    vals <- x[!is.na(x)]
    length(unique(vals)) > 1
  })
  if (any(!col_vars)) {
    n_mono <- sum(!col_vars)
    if (verbose) message(sprintf("  Removing %d monomorphic SNPs", n_mono))
    X_sub <- X_sub[, col_vars, drop = FALSE]
    aligned <- aligned[col_vars, , drop = FALSE]
    snp_ids_for_geno <- aligned$bim_id
  }

  n_missing <- sum(is.na(X_sub))
  if (verbose) message(sprintf("Computing LD for %d SNPs from %d samples (%d missing genotypes)...",
                                ncol(X_sub), nrow(X_sub), n_missing))

  # Compute LD matrix.
  # When match_original = TRUE, use the C++ GCTA-style computation that matches
  # the binary's bfileOperations.cpp exactly: per-pair missing data handling,
  # sequential sample accumulation, same floating-point operation order.
  # Otherwise use R's standard sample correlation.
  if (match_original) {
    N_kept <- (nrow(X_sub) %/% 4L) * 4L
    if (N_kept < nrow(X_sub)) {
      X_sub <- X_sub[seq_len(N_kept), , drop = FALSE]
    }
    LD_mat <- compute_LD_gcta_cpp(X_sub, ncpus = 1L)
    colnames(LD_mat) <- rownames(LD_mat) <- colnames(X_sub)
  } else {
    LD_mat <- compute_LD(X_sub, method = "sample")
  }

  # 7. Prepare sum_stat for dentist() — needs pos and z columns
  dentist_input <- data.frame(
    SNP = aligned$SNP,
    pos = aligned$pos,
    z = aligned$z,
    stringsAsFactors = FALSE
  )

  # 8. Run dentist
  if (verbose) message("Running R DENTIST implementation...")
  dentist_result <- dentist(
    sum_stat = dentist_input,
    LD_mat = LD_mat,
    nSample = nSample,
    window_size = window_size,
    pValueThreshold = pValueThreshold,
    propSVD = propSVD,
    gcControl = gcControl,
    nIter = nIter,
    gPvalueThreshold = gPvalueThreshold,
    duprThreshold = duprThreshold,
    use_original_windowing = use_original_windowing,
    min_dim = min_dim,
    ncpus = ncpus,
    seed = seed,
    correct_chen_et_al_bug = correct_chen_et_al_bug,
    match_original = match_original
  )

  # Attach SNP names to result using index_global when available (windowed mode),
  # or directly when the result has the same number of rows as the input.
  if ("index_global" %in% colnames(dentist_result)) {
    dentist_result$SNP <- dentist_input$SNP[dentist_result$index_global]
  } else if (nrow(dentist_result) == nrow(dentist_input)) {
    dentist_result$SNP <- dentist_input$SNP
  } else {
    warning("Cannot assign SNP names: row count mismatch and no index_global column")
  }

  if (verbose) {
    n_outlier <- sum(dentist_result$outlier, na.rm = TRUE)
    message(sprintf("Done: %d variants tested, %d outliers detected (%.1f%%)",
                    nrow(dentist_result), n_outlier,
                    100 * n_outlier / nrow(dentist_result)))
  }

  return(list(
    result = dentist_result,
    sum_stat = aligned,
    LD_mat = LD_mat
  ))
}

#' Parse DENTIST Binary Output Files
#'
#' Reads the output files produced by the DENTIST C++ binary and returns
#' a structured data frame. Useful for comparing binary output against
#' \code{\link{dentist_from_files}} results.
#'
#' @param output_prefix The output prefix used when running the DENTIST binary
#'   (the \code{--out} argument).
#' @param pValueThreshold P-value threshold used to determine outlier status.
#'   Default is 5e-8.
#' @return A data frame with columns: SNP, test_stat, neg_log10_p, is_duplicate, outlier.
#'
#' @details
#' The DENTIST binary (non-debug mode) writes a \code{.DENTIST.full.txt} file with
#' 4 tab-separated columns (no header): rsID, stat/lambda, -log10(pvalue), isDuplicate.
#' It also writes a \code{.DENTIST.short.txt} file listing outlier SNP IDs (one per line).
#'
#' @export
parse_dentist_output <- function(output_prefix, pValueThreshold = 5e-8) {
  full_file <- paste0(output_prefix, ".DENTIST.full.txt")
  short_file <- paste0(output_prefix, ".DENTIST.short.txt")

  if (!file.exists(full_file)) {
    stop(paste0("DENTIST full output file not found: ", full_file))
  }

  full_df <- read.table(full_file, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE,
                         col.names = c("SNP", "test_stat", "neg_log10_p", "is_duplicate"))
  full_df$is_duplicate <- as.logical(full_df$is_duplicate)
  full_df$outlier <- full_df$neg_log10_p > -log10(pValueThreshold)

  # Cross-check with short file if it exists
  if (file.exists(short_file) && file.info(short_file)$size > 0) {
    short_snps <- trimws(readLines(short_file))
    short_snps <- short_snps[nchar(short_snps) > 0]
    n_outlier_full <- sum(full_df$outlier)
    n_outlier_short <- length(short_snps)
    if (n_outlier_full != n_outlier_short) {
      warning(paste0("Outlier count mismatch: full file has ", n_outlier_full,
                      " outliers, short file has ", n_outlier_short, " outliers."))
    }
  }

  return(full_df)
}
