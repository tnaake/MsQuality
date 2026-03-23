## create test sets for test_function_Chromatograms_metrics.R
## create toy example chromatographic data using ChromBackendMemory

library(Chromatograms)

## Standard test chromatograms with 3 chromatograms (one empty)
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

## Chromatograms with 10x intensity changes for msSignal10xChange tests
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

## ============================================================
## Per-chromatogram scalar metrics
## ============================================================

test_that("chromatographyDuration returns per-chromatogram durations.", {
    tmp <- chromatographyDuration(chr)
    ## 3 chromatograms: durations are (3.9-2.1), NA (empty), (7.5-5.1)
    expect_equal(length(tmp), 3)
    expect_equal(as.numeric(tmp[1]), 1.8, tolerance = 1e-6)
    expect_true(is.na(tmp[2]))
    expect_equal(as.numeric(tmp[3]), 2.4, tolerance = 1e-6)
    expect_equal(attr(tmp, "chromatographyDuration"), "MS:4000053")
})

test_that("peakCount returns per-chromatogram data point counts.", {
    tmp <- peakCount(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], 5L)
    expect_equal(tmp[2], 0L)
    expect_equal(tmp[3], 5L)
})

test_that("numberEmptyScans works on Chromatograms.", {
    tmp <- numberEmptyScans(chr, msLevel = 1L)
    expect_equal(as.numeric(tmp), 1)
})

test_that("rtAcquisitionRange returns per-chromatogram matrix.", {
    tmp <- rtAcquisitionRange(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_equal(ncol(tmp), 2)
    expect_equal(colnames(tmp), c("min", "max"))
    expect_equal(tmp[1, ], c(min = 2.1, max = 3.9))
    expect_true(all(is.na(tmp[2, ])))
    expect_equal(tmp[3, ], c(min = 5.1, max = 7.5))
    expect_equal(attr(tmp, "rtAcquisitionRange"), "MS:4000070")
})

test_that("maxIntensity returns per-chromatogram values.", {
    tmp <- maxIntensity(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], 400)
    expect_true(is.na(tmp[2]))
    expect_equal(tmp[3], 1200)
})

test_that("intensityMean returns per-chromatogram values.", {
    tmp <- intensityMean(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], mean(c(100, 250, 400, 300, 150)))
    expect_true(is.na(tmp[2]))
    expect_equal(tmp[3], mean(c(80, 500, 1200, 600, 120)))
})

test_that("intensitySd returns per-chromatogram values.", {
    tmp <- intensitySd(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], sd(c(100, 250, 400, 300, 150)))
    expect_true(is.na(tmp[2]))
    expect_equal(tmp[3], sd(c(80, 500, 1200, 600, 120)))
})

test_that("intensityQuartiles returns per-chromatogram matrix.", {
    tmp <- intensityQuartiles(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_equal(ncol(tmp), 3)
    expect_equal(colnames(tmp), c("Q1", "Q2", "Q3"))
    ## Check chromatogram 1
    q1 <- unname(quantile(c(100, 250, 400, 300, 150), probs = c(0.25, 0.5, 0.75)))
    expect_equal(as.numeric(tmp[1, ]), q1, tolerance = 1e-6)
    ## Empty chromatogram -> all NA
    expect_true(all(is.na(tmp[2, ])))
    ## Chromatogram 3
    q3 <- unname(quantile(c(80, 500, 1200, 600, 120), probs = c(0.25, 0.5, 0.75)))
    expect_equal(as.numeric(tmp[3, ]), q3, tolerance = 1e-6)
})

test_that("intensityRange returns per-chromatogram matrix.", {
    tmp <- intensityRange(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_equal(ncol(tmp), 2)
    expect_equal(colnames(tmp), c("min", "max"))
    expect_equal(tmp[1, ], c(min = 100, max = 400))
    expect_true(all(is.na(tmp[2, ])))
    expect_equal(tmp[3, ], c(min = 80, max = 1200))
})

test_that("rtIqr returns per-chromatogram values.", {
    tmp <- rtIqr(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], IQR(c(2.1, 2.5, 3.0, 3.4, 3.9)), tolerance = 1e-6)
    expect_true(is.na(tmp[2]))
    expect_equal(tmp[3], IQR(c(5.1, 5.8, 6.3, 6.9, 7.5)), tolerance = 1e-6)
})

test_that("baselineIntensity returns per-chromatogram values.", {
    tmp <- baselineIntensity(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], as.numeric(quantile(c(100, 250, 400, 300, 150), 0.05)))
    expect_true(is.na(tmp[2]))
    expect_equal(tmp[3], as.numeric(quantile(c(80, 500, 1200, 600, 120), 0.05)))
})

test_that("signalToNoiseRatio returns per-chromatogram values.", {
    tmp <- signalToNoiseRatio(chr)
    expect_equal(length(tmp), 3)
    ## chromatogram 1 and 3 should have finite S/N > 0
    expect_true(is.finite(tmp[1]) && tmp[1] > 0)
    expect_true(is.na(tmp[2]))
    expect_true(is.finite(tmp[3]) && tmp[3] > 0)
})

test_that("msSignal10xChange returns per-chromatogram values on Chromatograms.", {
    tmp_jump <- msSignal10xChange(chr_jump, change = "jump", msLevel = 1L)
    tmp_fall <- msSignal10xChange(chr_jump, change = "fall", msLevel = 1L)
    expect_equal(length(tmp_jump), 2)
    expect_equal(length(tmp_fall), 2)
    ## chr_jump[1]: 50,100,2000,100,50 -> jump: 100->2000 = 20x (1 jump)
    ## chr_jump[2]: 200,20,30,20,200 -> jump: 20->200 = 10x (1 jump)
    expect_equal(tmp_jump[1], 1L)
    expect_equal(tmp_jump[2], 1L)
    ## chr_jump[1]: 2000->100 = 20x fall (1 fall)
    ## chr_jump[2]: 200->20 = 10x fall (1 fall)
    expect_equal(tmp_fall[1], 1L)
    expect_equal(tmp_fall[2], 1L)
})

test_that("medianTicRtIqr returns per-chromatogram values on Chromatograms.", {
    tmp <- medianTicRtIqr(chr, msLevel = 1L)
    expect_equal(length(tmp), 3)
    ## chromatogram 1: RTs=2.1,2.5,3.0,3.4,3.9; Q1=2.3, Q3=3.65
    rt1 <- c(2.1, 2.5, 3.0, 3.4, 3.9)
    int1 <- c(100, 250, 400, 300, 150)
    q1 <- quantile(rt1, 0.25)
    q3 <- quantile(rt1, 0.75)
    in_iqr <- int1[rt1 >= q1 & rt1 <= q3]
    expect_equal(tmp[1], median(in_iqr), tolerance = 1e-6)
    ## empty -> NA
    expect_true(is.na(tmp[2]))
    ## chromatogram 3
    rt3 <- c(5.1, 5.8, 6.3, 6.9, 7.5)
    int3 <- c(80, 500, 1200, 600, 120)
    q1_3 <- quantile(rt3, 0.25)
    q3_3 <- quantile(rt3, 0.75)
    in_iqr_3 <- int3[rt3 >= q1_3 & rt3 <= q3_3]
    expect_equal(tmp[3], median(in_iqr_3), tolerance = 1e-6)
})

## ============================================================
## Per-chromatogram EIC metrics (Part 2)
## ============================================================

test_that("xicFwhm returns per-chromatogram values.", {
    tmp <- xicFwhm(chr)
    expect_equal(length(tmp), 3)
    expect_true(tmp[1] > 0)
    expect_true(is.na(tmp[2]))  # empty
    expect_true(tmp[3] > 0)
})

test_that("xicFwhm handles NA intensities without error.", {
    cdata_na <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_na <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        intensity = c(NA, 20, 50, 100, 500, 100, 50, NA, 10, 5)
    ))
    chr_na <- Chromatograms(ChromBackendMemory(), chromData = cdata_na,
                            peaksData = pdata_na)
    expect_no_error(xicFwhm(chr_na))
    tmp <- xicFwhm(chr_na)
    expect_equal(length(tmp), 1)
})

test_that("xicFwhm accepts pre-computed peakBoundary matrix.", {
    pb <- peakBoundary(chr)
    tmp <- xicFwhm(chr, peakBoundary = pb)
    expect_equal(length(tmp), 3)
})

test_that("peakBoundary length mismatch raises error.", {
    ## Create a 1-row peakBoundary matrix (chr has 3 chromatograms)
    bad_pb <- matrix(c(1, 5), nrow = 1,
                     dimnames = list(NULL, c("left_boundary", "right_boundary")))
    expect_error(xicFwhm(chr, peakBoundary = bad_pb), "peakBoundary")
    expect_error(peakWidth(chr, peakBoundary = bad_pb), "peakBoundary")
    expect_error(gaussianSimilarity(chr, peakBoundary = bad_pb), "peakBoundary")
    expect_error(peakProminence(chr, peakBoundary = bad_pb), "peakBoundary")
})



test_that("peakWidth returns per-chromatogram values.", {
    tmp <- peakWidth(chr)
    expect_equal(length(tmp), 3)
    expect_true(tmp[1] > 0)
    expect_true(is.na(tmp[2]))
    expect_true(tmp[3] > 0)
})

test_that("peakWidth returns correct width for clean peak.", {
    cdata_pw <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pw <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_pw <- Chromatograms(ChromBackendMemory(), chromData = cdata_pw,
                            peaksData = pdata_pw)
    tmp <- peakWidth(chr_pw)
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 6)
})

test_that("peakWidth accepts pre-computed peakBoundary matrix.", {
    pb <- peakBoundary(chr)
    tmp <- peakWidth(chr, peakBoundary = pb)
    expect_equal(length(tmp), 3)
    expect_true(is.na(tmp[2]))
})

test_that("gaussianSimilarity returns per-chromatogram matrix.", {
    ## Create a multi-chrom object with enough points per chrom
    cdata_beta <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100.0, 200.0),
        dataOrigin = c("mem1", "mem1")
    )
    pdata_beta <- list(
        data.frame(
            rtime = seq(1, 20, by = 0.5),
            intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
                          700, 650, 500, 350, 200, 100, 50, 25, 15,
                          10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1)
        ),
        data.frame(
            rtime = c(1, 2, 3, 4, 5, 6, 7),
            intensity = c(0, 10, 50, 100, 50, 10, 0)
        )
    )
    chr_beta <- Chromatograms(ChromBackendMemory(), chromData = cdata_beta,
                              peaksData = pdata_beta)
    tmp <- gaussianSimilarity(chr_beta)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 2)
    expect_equal(ncol(tmp), 2)
    expect_equal(colnames(tmp), c("gaussian_similarity", "gaussian_residuals"))
    ## First chromatogram (39 points, clear peak) should have valid values
    expect_true(!is.na(tmp[1, "gaussian_similarity"]))
    ## Second chromatogram (7 points) may or may not have enough points
})

test_that("gaussianSimilarity returns NA for empty chromatogram.", {
    cdata_empty <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_empty <- list(data.frame(rtime = numeric(), intensity = numeric()))
    chr_empty <- Chromatograms(ChromBackendMemory(), chromData = cdata_empty,
                               peaksData = pdata_empty)
    tmp <- gaussianSimilarity(chr_empty)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 1)
    expect_true(all(is.na(tmp)))
})

test_that("gaussianSimilarity returns NA for peak with < 5 points.", {
    cdata_short <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_short <- list(data.frame(
        rtime = c(1, 2, 3),
        intensity = c(10, 100, 10)
    ))
    chr_short <- Chromatograms(ChromBackendMemory(), chromData = cdata_short,
                               peaksData = pdata_short)
    tmp <- gaussianSimilarity(chr_short)
    expect_true(all(is.na(tmp)))
})

test_that("gaussianSimilarity accepts pre-computed peakBoundary matrix.", {
    cdata_beta <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_beta <- list(data.frame(
        rtime = seq(1, 20, by = 0.5),
        intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
                      700, 650, 500, 350, 200, 100, 50, 25, 15,
                      10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1)
    ))
    chr_beta <- Chromatograms(ChromBackendMemory(), chromData = cdata_beta,
                              peaksData = pdata_beta)
    pb <- peakBoundary(chr_beta)
    tmp <- gaussianSimilarity(chr_beta, peakBoundary = pb)
    expect_true(is.matrix(tmp))
    expect_true(!is.na(tmp[1, "gaussian_similarity"]))
})

test_that("peakProminence returns per-chromatogram values.", {
    tmp <- peakProminence(chr)
    expect_equal(length(tmp), 3)
    expect_true(tmp[1] > 0)
    expect_true(is.na(tmp[2]))  # empty
    expect_true(tmp[3] > 0)
})

test_that("peakProminence returns NA for zero-baseline chromatogram.", {
    cdata_zero <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_zero <- list(data.frame(
        rtime = seq(1, 10, by = 1),
        intensity = c(0, 0, 0, 0, 100, 0, 0, 0, 0, 0)
    ))
    chr_zero <- Chromatograms(ChromBackendMemory(), chromData = cdata_zero,
                              peaksData = pdata_zero)
    tmp <- peakProminence(chr_zero)
    expect_true(is.na(tmp))
})

test_that("peakProminence accepts pre-computed peakBoundary matrix.", {
    pb <- peakBoundary(chr)
    tmp <- peakProminence(chr, peakBoundary = pb)
    expect_equal(length(tmp), 3)
})

## ============================================================
## Shared metrics (from function_Spectra_metrics.R)
## ============================================================

test_that("ticQuantileRtFraction returns per-chromatogram matrix.", {
    tmp <- ticQuantileRtFraction(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_equal(ncol(tmp), 5)
    expect_equal(colnames(tmp), c("0%", "25%", "50%", "75%", "100%"))
    ## Row 1: valid fractions (0 to 1)
    expect_equal(unname(tmp[1, "0%"]), 0, tolerance = 1e-6)
    expect_equal(unname(tmp[1, "100%"]), 1, tolerance = 1e-6)
    ## Row 2: empty → all NA
    expect_true(all(is.na(tmp[2, ])))
    ## Row 3: valid fractions
    expect_equal(unname(tmp[3, "0%"]), 0, tolerance = 1e-6)
    expect_equal(unname(tmp[3, "100%"]), 1, tolerance = 1e-6)
    expect_equal(attr(tmp, "ticQuantileRtFraction"), "MS:4000183")
})

test_that("areaUnderTic returns per-chromatogram values.", {
    tmp <- areaUnderTic(chr)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], sum(c(100, 250, 400, 300, 150)))
    expect_equal(tmp[2], 0)  # empty chromatogram: sum of nothing = 0
    expect_equal(tmp[3], sum(c(80, 500, 1200, 600, 120)))
    expect_equal(attr(tmp, "areaUnderTic"), "MS:4000155")
})

test_that("areaUnderTicRtQuantiles returns per-chromatogram matrix.", {
    tmp <- areaUnderTicRtQuantiles(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_equal(ncol(tmp), 4)
    expect_equal(colnames(tmp), c("25%", "50%", "75%", "100%"))
    ## Row 1: valid positive areas
    expect_true(all(tmp[1, ] >= 0))
    ## Row 2: empty → NA
    expect_true(all(is.na(tmp[2, ])))
    ## Row 3: valid positive areas
    expect_true(all(tmp[3, ] >= 0))
    expect_equal(attr(tmp, "areaUnderTicRtQuantiles"), "MS:4000156")
})

## Chromatograms with MS1 and MS2 for msLevel filter tests
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

test_that("areaUnderTic filters by msLevel per chromatogram.", {
    tmp_ms1 <- areaUnderTic(chr_ms, msLevel = 1L)
    ## MS1 chromatograms: chr1=600, chr3=1500
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1[1], 600)
    expect_equal(tmp_ms1[2], 1500)
    expect_equal(attr(tmp_ms1, "areaUnderTic"), "MS:4000155")

    tmp_ms2 <- areaUnderTic(chr_ms, msLevel = 2L)
    ## MS2 chromatograms: chr2=180
    expect_equal(length(tmp_ms2), 1)
    expect_equal(tmp_ms2[1], 180)
    expect_equal(attr(tmp_ms2, "areaUnderTic"), "MS:4000155")
})

test_that("rtAcquisitionRange filters by msLevel.", {
    tmp_ms1 <- rtAcquisitionRange(chr_ms, msLevel = 1L)
    expect_true(is.matrix(tmp_ms1))
    expect_equal(nrow(tmp_ms1), 2)  # 2 MS1 chromatograms
})

test_that("rtIqr filters by msLevel.", {
    tmp_ms1 <- rtIqr(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)  # 2 MS1 chromatograms
})



test_that("numberEmptyScans filters by msLevel.", {
    tmp <- numberEmptyScans(chr_ms, msLevel = 1L)
    expect_equal(as.numeric(tmp), 0)
})

test_that("msSignal10xChange filters by msLevel.", {
    tmp <- msSignal10xChange(chr_ms, change = "jump", msLevel = 1L)
    expect_equal(length(tmp), 2)  # 2 MS1 chromatograms
})

## ============================================================
## Metrics with NA intensities
## ============================================================

test_that("scalar metrics handle NA values correctly per-chromatogram.", {
    cdata_na <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(112.2, 123.3),
        dataOrigin = c("mem1", "mem1")
    )
    pdata_na <- list(
        data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
                   intensity = c(100, NA, 400, NA, 150)),
        data.frame(rtime = c(5.1, 5.8, 6.3),
                   intensity = c(80, 500, NA))
    )
    chr_na <- Chromatograms(ChromBackendMemory(), chromData = cdata_na,
                            peaksData = pdata_na)

    tmp_max <- maxIntensity(chr_na, na.rm = TRUE)
    expect_equal(tmp_max[1], 400)
    expect_equal(tmp_max[2], 500)

    tmp_mean <- intensityMean(chr_na, na.rm = TRUE)
    expect_equal(tmp_mean[1], mean(c(100, 400, 150)))
    expect_equal(tmp_mean[2], mean(c(80, 500)))

    tmp_range <- intensityRange(chr_na, na.rm = TRUE)
    expect_equal(tmp_range[1, ], c(min = 100, max = 400))
    expect_equal(tmp_range[2, ], c(min = 80, max = 500))
})

## ============================================================
## All-zero chromatograms
## ============================================================

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

test_that("maxIntensity handles all-zero chromatograms.", {
    tmp <- maxIntensity(chr_zero)
    expect_equal(length(tmp), 2)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

test_that("intensityMean handles all-zero chromatograms.", {
    tmp <- intensityMean(chr_zero)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

test_that("intensitySd handles all-zero chromatograms.", {
    tmp <- intensitySd(chr_zero)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

test_that("intensityQuartiles handles all-zero chromatograms.", {
    tmp <- intensityQuartiles(chr_zero)
    expect_true(all(tmp == 0))
})

test_that("intensityRange handles all-zero chromatograms.", {
    tmp <- intensityRange(chr_zero)
    expect_equal(tmp[1, ], c(min = 0, max = 0))
    expect_equal(tmp[2, ], c(min = 0, max = 0))
})

test_that("peakCount handles all-zero chromatograms.", {
    tmp <- peakCount(chr_zero)
    expect_equal(tmp[1], 5L)
    expect_equal(tmp[2], 7L)
})

test_that("baselineIntensity handles all-zero chromatograms.", {
    tmp <- baselineIntensity(chr_zero)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

test_that("signalToNoiseRatio handles all-zero chromatograms.", {
    tmp <- signalToNoiseRatio(chr_zero)
    ## MAD of all zeros = 0, so should return NA
    expect_true(is.na(tmp[1]))
    expect_true(is.na(tmp[2]))
})

test_that("xicFwhm handles all-zero chromatograms.", {
    tmp <- xicFwhm(chr_zero)
    expect_true(is.na(tmp[1]))
    expect_true(is.na(tmp[2]))
})

test_that("peakWidth handles all-zero chromatograms.", {
    tmp <- peakWidth(chr_zero)
    expect_true(is.na(tmp[1]))
    expect_true(is.na(tmp[2]))
})

test_that("gaussianSimilarity handles all-zero chromatograms.", {
    tmp <- gaussianSimilarity(chr_zero)
    expect_true(all(is.na(tmp)))
})

test_that("peakProminence handles all-zero chromatograms.", {
    tmp <- peakProminence(chr_zero)
    expect_true(is.na(tmp[1]))
    expect_true(is.na(tmp[2]))
})

test_that("msSignal10xChange handles all-zero chromatograms.", {
    tmp <- msSignal10xChange(chr_zero, change = "jump", msLevel = 1L)
    expect_equal(tmp[1], 0L)
    expect_equal(tmp[2], 0L)
})

test_that("areaUnderTic handles all-zero chromatograms.", {
    tmp <- areaUnderTic(chr_zero)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

test_that("medianTicRtIqr handles all-zero chromatograms.", {
    tmp <- medianTicRtIqr(chr_zero, msLevel = 1L)
    expect_equal(tmp[1], 0)
    expect_equal(tmp[2], 0)
})

## ============================================================
## Chromatograms with NA values
## ============================================================

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

test_that("maxIntensity handles NA chromatograms.", {
    tmp <- maxIntensity(chr_na2)
    expect_equal(tmp[1], 400)
    expect_true(is.na(tmp[2]))  # all-NA
    expect_equal(tmp[3], 100)
})

test_that("intensityMean handles NA chromatograms.", {
    tmp <- intensityMean(chr_na2)
    expect_equal(tmp[1], mean(c(100, 400, 150)))
    expect_true(is.nan(tmp[2]))  # mean of all-NA with na.rm = TRUE
    expect_equal(tmp[3], mean(c(10, 50, 100, 50, 10)))
})

test_that("intensitySd handles NA chromatograms.", {
    tmp <- intensitySd(chr_na2)
    expect_equal(tmp[1], sd(c(100, 400, 150)))
    expect_true(is.na(tmp[2]))  # sd of all-NA
    expect_equal(tmp[3], sd(c(10, 50, 100, 50, 10)))
})

test_that("intensityQuartiles handles NA chromatograms.", {
    tmp <- intensityQuartiles(chr_na2)
    expect_equal(nrow(tmp), 3)
    ## chromatograms with non-NA values should have valid quartiles
    expect_true(!is.na(tmp[1, "Q3"]))
    expect_true(!is.na(tmp[3, "Q3"]))
})

test_that("intensityRange handles NA chromatograms.", {
    tmp <- intensityRange(chr_na2)
    expect_equal(tmp[1, ], c(min = 100, max = 400))
    expect_true(all(is.na(tmp[2, ])))  # all-NA
    expect_equal(tmp[3, ], c(min = 10, max = 100))
})

test_that("peakCount handles NA chromatograms.", {
    tmp <- peakCount(chr_na2, na.rm = FALSE)
    expect_equal(tmp[1], 5L)  # counts NAs
    expect_equal(tmp[2], 5L)
    expect_equal(tmp[3], 7L)

    tmp_rm <- peakCount(chr_na2, na.rm = TRUE)
    expect_equal(tmp_rm[1], 3L)  # only non-NA
    expect_equal(tmp_rm[2], 0L)
    expect_equal(tmp_rm[3], 5L)
})

test_that("baselineIntensity handles NA chromatograms.", {
    tmp <- baselineIntensity(chr_na2)
    expect_equal(tmp[1], as.numeric(quantile(c(100, NA, 400, NA, 150),
                                             0.05, na.rm = TRUE)))
    expect_true(!is.na(tmp[3]))
})

test_that("signalToNoiseRatio handles NA chromatograms.", {
    tmp <- signalToNoiseRatio(chr_na2)
    expect_true(is.na(tmp[2]))  # all-NA
    ## chromatograms 1 and 3 should produce a value or NA, but not error
    expect_equal(length(tmp), 3)
})

test_that("xicFwhm handles NA chromatograms.", {
    tmp <- xicFwhm(chr_na2)
    expect_equal(length(tmp), 3)
    expect_true(is.na(tmp[2]))  # all-NA
})

test_that("peakWidth handles NA chromatograms.", {
    tmp <- peakWidth(chr_na2)
    expect_equal(length(tmp), 3)
    expect_true(is.na(tmp[2]))  # all-NA
})

test_that("gaussianSimilarity handles NA chromatograms.", {
    tmp <- gaussianSimilarity(chr_na2)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), 3)
    expect_true(all(is.na(tmp[2, ])))  # all-NA
})

test_that("peakProminence handles NA chromatograms.", {
    tmp <- peakProminence(chr_na2)
    expect_equal(length(tmp), 3)
    expect_true(is.na(tmp[2]))  # all-NA
})

test_that("msSignal10xChange handles NA chromatograms.", {
    tmp <- msSignal10xChange(chr_na2, change = "jump", msLevel = 1L)
    expect_equal(length(tmp), 3)
})

test_that("areaUnderTic handles NA chromatograms.", {
    tmp <- areaUnderTic(chr_na2)
    expect_equal(length(tmp), 3)
    expect_equal(tmp[1], sum(c(100, 400, 150)))  # na.rm = TRUE
})

test_that("medianTicRtIqr handles NA chromatograms.", {
    tmp <- medianTicRtIqr(chr_na2, msLevel = 1L)
    expect_equal(length(tmp), 3)
    expect_true(is.na(tmp[2]))  # all-NA
})

## ============================================================
## msLevel filtering for intensity-related metrics
## ============================================================

## Re-use chr_ms: 3 chromatograms with msLevel = c(1L, 2L, 1L)
## chr_ms pdata: chr1(100,200,300), chr2(50,60,70), chr3(400,500,600)

test_that("maxIntensity filters by msLevel.", {
    tmp_all <- maxIntensity(chr_ms)
    expect_equal(length(tmp_all), 3)
    tmp_ms1 <- maxIntensity(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1, c(300, 600))
    tmp_ms2 <- maxIntensity(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
    expect_equal(tmp_ms2, 70)
})

test_that("intensityMean filters by msLevel.", {
    tmp_ms1 <- intensityMean(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1[1], mean(c(100, 200, 300)))
    expect_equal(tmp_ms1[2], mean(c(400, 500, 600)))
    tmp_ms2 <- intensityMean(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
    expect_equal(tmp_ms2[1], mean(c(50, 60, 70)))
})

test_that("intensitySd filters by msLevel.", {
    tmp_ms1 <- intensitySd(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1[1], sd(c(100, 200, 300)))
    tmp_ms2 <- intensitySd(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
    expect_equal(tmp_ms2[1], sd(c(50, 60, 70)))
})

test_that("intensityQuartiles filters by msLevel.", {
    tmp_ms1 <- intensityQuartiles(chr_ms, msLevel = 1L)
    expect_true(is.matrix(tmp_ms1))
    expect_equal(nrow(tmp_ms1), 2)
    expect_equal(ncol(tmp_ms1), 3)
    q_chr1 <- unname(quantile(c(100, 200, 300), probs = c(0.25, 0.5, 0.75)))
    expect_equal(as.numeric(tmp_ms1[1, ]), q_chr1, tolerance = 1e-6)
    tmp_ms2 <- intensityQuartiles(chr_ms, msLevel = 2L)
    expect_equal(nrow(tmp_ms2), 1)
})

test_that("intensityRange filters by msLevel.", {
    tmp_ms1 <- intensityRange(chr_ms, msLevel = 1L)
    expect_true(is.matrix(tmp_ms1))
    expect_equal(nrow(tmp_ms1), 2)
    expect_equal(tmp_ms1[1, ], c(min = 100, max = 300))
    expect_equal(tmp_ms1[2, ], c(min = 400, max = 600))
    tmp_ms2 <- intensityRange(chr_ms, msLevel = 2L)
    expect_equal(nrow(tmp_ms2), 1)
    expect_equal(tmp_ms2[1, ], c(min = 50, max = 70))
})

test_that("peakCount filters by msLevel.", {
    tmp_ms1 <- peakCount(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1, c(3L, 3L))
    tmp_ms2 <- peakCount(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
    expect_equal(tmp_ms2, 3L)
})

test_that("baselineIntensity filters by msLevel.", {
    tmp_ms1 <- baselineIntensity(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    expect_equal(tmp_ms1[1],
                 as.numeric(quantile(c(100, 200, 300), 0.05)))
    tmp_ms2 <- baselineIntensity(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
})

test_that("signalToNoiseRatio filters by msLevel.", {
    tmp_ms1 <- signalToNoiseRatio(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    tmp_ms2 <- signalToNoiseRatio(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
})

test_that("xicFwhm filters by msLevel.", {
    tmp_ms1 <- xicFwhm(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    tmp_ms2 <- xicFwhm(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
})

test_that("peakWidth filters by msLevel.", {
    tmp_ms1 <- peakWidth(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    tmp_ms2 <- peakWidth(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
})

test_that("gaussianSimilarity filters by msLevel.", {
    tmp_ms1 <- gaussianSimilarity(chr_ms, msLevel = 1L)
    expect_true(is.matrix(tmp_ms1))
    expect_equal(nrow(tmp_ms1), 2)
    tmp_ms2 <- gaussianSimilarity(chr_ms, msLevel = 2L)
    expect_equal(nrow(tmp_ms2), 1)
})

test_that("peakProminence filters by msLevel.", {
    tmp_ms1 <- peakProminence(chr_ms, msLevel = 1L)
    expect_equal(length(tmp_ms1), 2)
    tmp_ms2 <- peakProminence(chr_ms, msLevel = 2L)
    expect_equal(length(tmp_ms2), 1)
})

test_that("msLevel filtering with empty result returns zero-length output.", {
    ## No msLevel 3 chromatograms exist
    tmp <- maxIntensity(chr_ms, msLevel = 3L)
    expect_equal(length(tmp), 0)
    tmp_q <- intensityQuartiles(chr_ms, msLevel = 3L)
    expect_equal(nrow(tmp_q), 0)
    tmp_r <- intensityRange(chr_ms, msLevel = 3L)
    expect_equal(nrow(tmp_r), 0)
})

test_that("msLevel = integer() (default) returns all chromatograms.", {
    tmp_default <- maxIntensity(chr_ms)
    tmp_explicit <- maxIntensity(chr_ms, msLevel = integer())
    expect_equal(tmp_default, tmp_explicit)
})

