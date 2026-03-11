# Extracted from test_function_Chromatograms_metrics.R:511

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
cdata_ms <- data.frame(
    msLevel = c(1L, 2L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem1", "mem1")
)
pdata_ms <- list(
    data.frame(rtime = c(2.1, 2.5, 3.0), intensity = c(100, 200, 300)),
    data.frame(rtime = c(3.5, 4.0, 4.5), intensity = c(50, 60, 70)),
    data.frame(rtime = c(5.0, 5.5, 6.0), intensity = c(400, 500, 600))
)
chr_ms <- Chromatograms(ChromBackendMemory(), chromData = cdata_ms,
                        peaksData = pdata_ms)

# test -------------------------------------------------------------------------
tmp_ms1 <- areaUnderTic(chr_ms, msLevel = 1L)
expect_equal(length(tmp_ms1), 2)
expect_equal(tmp_ms1[1], 600)
expect_equal(tmp_ms1[2], 1500)
expect_equal(attr(tmp_ms1, "areaUnderTic"), "MS:4000029")
tmp_ms2 <- areaUnderTic(chr_ms, msLevel = 2L)
expect_equal(length(tmp_ms2), 1)
expect_equal(tmp_ms2[1], 180)
expect_equal(attr(tmp_ms2, "areaUnderTic"), "MS:4000030")
