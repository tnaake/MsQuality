# Extracted from test_function_Chromatograms_metrics.R:467

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MsQuality", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(Chromatograms)
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem1", "mem1")
)
pdata <- list(
    data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
               intensity = c(100, 250, 400, 300, 150)),
    data.frame(rtime = numeric(), intensity = numeric()),
    data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
               intensity = c(80, 500, 1200, 600, 120))
)
chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
cdata_jump <- data.frame(
    msLevel = c(1L, 1L),
    mz = c(112.2, 123.3),
    dataOrigin = c("mem1", "mem1")
)
pdata_jump <- list(
    data.frame(rtime = c(1.0, 2.0, 3.0, 4.0, 5.0),
               intensity = c(50, 100, 2000, 100, 50)),
    data.frame(rtime = c(1.0, 2.0, 3.0, 4.0, 5.0),
               intensity = c(200, 20, 30, 20, 200))
)
chr_jump <- Chromatograms(ChromBackendMemory(), chromData = cdata_jump,
                          peaksData = pdata_jump)

# test -------------------------------------------------------------------------
tmp <- areaUnderTic(chr)
expect_equal(length(tmp), 3)
expect_equal(tmp[1], sum(c(100, 250, 400, 300, 150)))
expect_equal(tmp[2], 0)
expect_equal(tmp[3], sum(c(80, 500, 1200, 600, 120)))
expect_equal(attr(tmp, "areaUnderTic"), "MS:4000029")
