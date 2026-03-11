# Extracted from test_function_Spectra_metrics.R:258

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MsQuality", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
sciex_file <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
    MsDataHub::X20171016_POOL_POS_3_105.134.mzML())
sps_sciex <- Spectra(sciex_file)
set.seed(1)
sps_sciex@backend$precursorCharge <- as.integer(
    sample(x = c(1, 2, 3, 4), size = 1862,
           replace = TRUE, prob = c(0.7, 0.15, 0.1, 0.05)))
sps_sciex@backend$precursorMz <- rnorm(n = 1862, mean = 500, sd = 100)
sps_sciex@backend$precursorIntensity <- rpois(n = 1862, lambda = 10000)

# test -------------------------------------------------------------------------
expect_error(numberEmptyScans(NULL))
expect_error(numberEmptyScans(1:10))
expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 1L)), 0)
expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 2L)), 0)
expect_equal(as.numeric(numberEmptyScans(sps_sciex, msLevel = 3L)), 0)
spd <- DataFrame(
        msLevel = c(2L, 2L), polarity = c(1L, 1L),
        id = c("unknown", "HMDB0000001"),
        name = c("unknown", "1-Methylhistidine"))
