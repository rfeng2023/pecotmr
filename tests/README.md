# Testing pecotmr

## Running Tests

```r
# Run all tests
devtools::test()

# Run a specific test file
devtools::load_all()
testthat::test_file("tests/testthat/test_misc.R")
```

**Note:** Always use `devtools::test()` (not `testthat::test_dir()`) because it
properly loads the package namespace before running tests.

## Measuring Code Coverage (R only, excluding C++)

The package includes C++ code (via Rcpp/RcppArmadillo) that causes `covr` to
fail when `gcov` cannot parse compiled objects. To measure R-only coverage,
create a fake `gcov` wrapper that skips C++ instrumentation:

```bash
# 1. Create a fake gcov that exits silently
mkdir -p /tmp/fake_gcov
echo '#!/bin/bash' > /tmp/fake_gcov/gcov
echo 'exit 0' >> /tmp/fake_gcov/gcov
chmod +x /tmp/fake_gcov/gcov

# 2. Run coverage with the fake gcov on PATH
PATH="/tmp/fake_gcov:$PATH" Rscript -e '
cov <- covr::package_coverage(type = "tests", quiet = TRUE)

# Per-file coverage (proportion of expressions hit)
df <- as.data.frame(cov)
file_cov <- aggregate(value ~ filename, data = df,
                      FUN = function(x) round(100 * mean(x > 0), 1))
file_cov <- file_cov[order(file_cov$value), ]
for (i in seq_len(nrow(file_cov))) {
  cat(sprintf("%-45s %5.1f%%\n", file_cov$filename[i], file_cov$value[i]))
}
cat(sprintf("\nOverall R coverage: %.1f%%\n",
            round(100 * mean(df$value > 0), 1)))
'
```

**Key detail:** The `value` column in `covr`'s data frame is an execution count
(not a boolean). Use `mean(x > 0)` to compute the proportion of covered
expressions, not `mean(x)`.
