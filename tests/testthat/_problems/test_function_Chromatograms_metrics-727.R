# Extracted from test_function_Chromatograms_metrics.R:727

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
cdata_zero <- data.frame(
    msLevel = c(1L, 1L),
    mz = c(112.2, 123.3),
    dataOrigin = c("mem1", "mem1")
)
pdata_zero <- list(
    data.frame(rtime = c(1, 2, 3, 4, 5),
               intensity = c(0, 0, 0, 0, 0)),
    data.frame(rtime = c(6, 7, 8, 9, 10, 11, 12),
               intensity = c(0, 0, 0, 0, 0, 0, 0))
)
chr_zero <- Chromatograms(ChromBackendMemory(), chromData = cdata_zero,
                          peaksData = pdata_zero)
cdata_na2 <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem1", "mem1")
)
pdata_na2 <- list(
    ## Some NAs interspersed
    data.frame(rtime = c(1, 2, 3, 4, 5),
               intensity = c(100, NA, 400, NA, 150)),
    ## All NAs
    data.frame(rtime = c(6, 7, 8, 9, 10),
               intensity = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)),
    ## Leading/trailing NAs
    data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
               intensity = c(NA, 10, 50, 100, 50, 10, NA))
)
chr_na2 <- Chromatograms(ChromBackendMemory(), chromData = cdata_na2,
                         peaksData = pdata_na2)

# test -------------------------------------------------------------------------
tmp <- intensityQuartiles(chr_na2)
