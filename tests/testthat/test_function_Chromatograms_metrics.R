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

## Chromatograms with 10x intensity changes for intensity10xChange tests
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

test_that("chromatographyDuration works properly for Chromatograms.", {
    tmp <- chromatographyDuration(chr)
    ## Returns aggregated duration across all chromatograms: max(7.5) - min(2.1) = 5.4
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 5.4, tolerance = 1e-6)

    ## test attributes
    expect_equal(attr(tmp, "chromatographyDuration"), "MS:4000053")
})

test_that("chromatogramCount works properly.", {
    tmp <- chromatogramCount(chr)
    ## Returns number of chromatograms (3 in test data)
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 3)

    ## test attributes
    expect_equal(attr(tmp, "chromatogramCount"), "MS:4000071")
})

test_that("numberEmptyChrom works properly.", {
    tmp <- numberEmptyChrom(chr)
    expect_equal(as.numeric(tmp), 1)  # One empty chromatogram
})

test_that("rtAcquisitionRange works properly for Chromatograms.", {
    tmp <- rtAcquisitionRange(chr)
    expect_equal(as.numeric(tmp), c(2.1, 7.5))
    expect_equal(names(tmp), c("min", "max"))

    ## test attributes
    expect_equal(attr(tmp, "rtAcquisitionRange"), "MS:4000070")
})

test_that("maxIntensity works properly.", {
    tmp <- maxIntensity(chr)
    ## Aggregated max across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 1200)

    ## test attributes
    expect_equal(attr(tmp, "maxIntensity"), "custom_metric:max_intensity")
})

test_that("intensityMean works properly.", {
    tmp <- intensityMean(chr)
    ## Aggregated mean across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), mean(c(100, 250, 400, 300, 150, 80, 500, 1200, 600, 120)))

    ## test attributes
    expect_equal(attr(tmp, "intensityMean"), "custom_metric:intensity_mean")
})

test_that("intensitySd works properly.", {
    tmp <- intensitySd(chr)
    ## Aggregated standard deviation across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), sd(c(100, 250, 400, 300, 150, 80, 500, 1200, 600, 120)))

    ## test attributes
    expect_equal(attr(tmp, "intensitySd"), "custom_metric:intensity_sd")
})

test_that("intensityQuartiles works properly.", {
    tmp <- intensityQuartiles(chr)
    expect_equal(length(tmp), 6)
    expect_equal(names(tmp), c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"))

    ## test attributes
    expect_equal(attr(tmp, "intensityQuartiles"), "custom_metric:intensity_quartiles")
})

test_that("intensityRange works properly.", {
    tmp <- intensityRange(chr)
    expect_equal(length(tmp), 2)
    expect_equal(names(tmp), c("min", "max"))
    expect_equal(as.numeric(tmp), c(80, 1200))

    ## test attributes
    expect_equal(attr(tmp, "intensityRange"), "custom_metric:intensity_range")
})

test_that("peakCount works properly.", {
    tmp <- peakCount(chr)
    expect_equal(as.numeric(tmp[1]), 5)
    expect_equal(as.numeric(tmp[2]), 0)  # Empty chromatogram
    expect_equal(as.numeric(tmp[3]), 5)

    ## test attributes
    expect_equal(attr(tmp, "peakCount"), "custom_metric:peak_count")
})

test_that("peakCount works with na.rm parameter.", {
    ## Create chromatograms with NA values in intensities
    cdata_na <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(112.2, 123.3),
        dataOrigin = c("mem1", "mem1")
    )
    pdata_na <- list(
        data.frame(rtime = c(1.0, 2.0, 3.0, 4.0, 5.0),
                   intensity = c(100, NA, 300, NA, 500)),
        data.frame(rtime = c(1.0, 2.0, 3.0),
                   intensity = c(NA, NA, NA))
    )
    chr_na <- Chromatograms(ChromBackendMemory(), chromData = cdata_na,
                            peaksData = pdata_na)

    ## Without na.rm (default = FALSE): counts all data points including NA
    tmp_with_na <- peakCount(chr_na, na.rm = FALSE)
    expect_equal(as.numeric(tmp_with_na[1]), 5)  # All 5 points counted
    expect_equal(as.numeric(tmp_with_na[2]), 3)  # All 3 points counted

    ## With na.rm = TRUE: counts only non-NA intensity values
    tmp_no_na <- peakCount(chr_na, na.rm = TRUE)
    ## na.rm = TRUE returns total count of non-NA values across all chromatograms
    expect_true(is.numeric(tmp_no_na))
    ## chr1 has 3 non-NA values (100, 300, 500), chr2 has 0 non-NA values
    expect_equal(as.numeric(tmp_no_na), 3)

    ## test attributes preserved
    expect_equal(attr(tmp_no_na, "peakCount"), "custom_metric:peak_count")
})

test_that("rtIqrChromatograms works properly.", {
    tmp <- rtIqrChromatograms(chr)
    ## Aggregated IQR across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), IQR(c(2.1, 2.5, 3.0, 3.4, 3.9, 5.1, 5.8, 6.3, 6.9, 7.5)))
})

test_that("baselineIntensity works properly.", {
    tmp <- baselineIntensity(chr)
    ## Aggregated baseline across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), as.numeric(quantile(c(100, 250, 400, 300, 150, 80, 500, 1200, 600, 120), probs = 0.05)))

    ## test attributes
    expect_equal(attr(tmp, "baselineIntensity"), "custom_metric:baseline_intensity")
})

test_that("signalToNoiseRatio works properly.", {
    tmp <- signalToNoiseRatio(chr)
    ## All intensities: 100, 250, 400, 300, 150, 80, 500, 1200, 600, 120
    ## Signal = max = 1200
    ## Noise = median of non-zero = median(80, 100, 120, 150, 250, 300, 400, 500, 600, 1200) = 275
    ## SNR = 1200 / 275 = 4.363636
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 1200 / 275, tolerance = 1e-6)

    ## test attributes
    expect_equal(attr(tmp, "signalToNoiseRatio"), "custom_metric:signal_to_noise_ratio")
})

test_that("intensity10xChange works properly.", {
    tmp_jump <- intensity10xChange(chr_jump, change = "jump")
    tmp_fall <- intensity10xChange(chr_jump, change = "fall")

    ## Aggregated count across all chromatograms
    expect_equal(length(tmp_jump), 1)
    expect_equal(length(tmp_fall), 1)
    expect_equal(as.numeric(tmp_jump), 2)  # 20x jump (2000/100) + 10x jump (200/20)
    expect_equal(as.numeric(tmp_fall), 2)  # 20x fall (100/2000=0.05) + 10x fall (20/200=0.1)

    ## test attributes
    expect_equal(attr(tmp_jump, "intensity10xJump"), "custom_metric:intensity_10x_jump")
    expect_equal(attr(tmp_fall, "intensity10xFall"), "custom_metric:intensity_10x_fall")
})

test_that("medianIntensityRtIqr works properly.", {
    tmp <- medianIntensityRtIqr(chr)
    ## All RTs: 2.1, 2.5, 3.0, 3.4, 3.9, 5.1, 5.8, 6.3, 6.9, 7.5
    ## R's default type=7: Q1 RT = 3.1, Q3 RT = 6.175
    ## Points within RT IQR: (3.4, 300), (3.9, 150), (5.1, 80), (5.8, 500)
    ## Median of intensities [80, 150, 300, 500] = (150 + 300) / 2 = 225
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 225, tolerance = 1e-6)

    ## test attributes
    expect_equal(attr(tmp, "medianIntensityRtIqr"), "custom_metric:median_intensity_rt_iqr")
})

test_that("areaUnderIntensityRtQuantiles works properly.", {
    tmp <- areaUnderIntensityRtQuantiles(chr)
    ## RT range [2.1, 7.5], 4 equal bins with cuts at 3.45, 4.8, 6.15
    ## Uses trapezoidal integration with linear interpolation at boundaries
    ## Just verify structure and that values are reasonable (positive areas)
    expect_equal(length(tmp), 4)
    expect_equal(names(tmp), c("Q1", "Q2", "Q3", "Q4"))
    expect_true(all(is.numeric(tmp)))
    expect_true(all(tmp >= 0))  # Areas should be non-negative
    expect_true(sum(tmp) > 0)   # Total area should be positive

    ## test attributes
    expect_equal(attr(tmp, "areaUnderIntensityRtQuantiles"), "custom_metric:area_under_intensity_rt_quantiles")
})

test_that("xicFwhm works properly.", {
    tmp <- xicFwhm(chr)
    ## Returns FWHM for each chromatogram
    ## chr1: has data with peak, should return numeric FWHM
    ## chr2: empty, should return NA
    ## chr3: has data with peak, should return numeric FWHM
    expect_equal(length(tmp), 3)
    expect_true(is.numeric(tmp[1]))
    expect_true(is.na(tmp[2]))  # Empty chromatogram
    expect_true(is.numeric(tmp[3]))
    expect_true(tmp[1] > 0)  # FWHM should be positive
    expect_true(tmp[3] > 0)  # FWHM should be positive

    ## test attributes
    expect_equal(attr(tmp, "xicFwhm"), "custom_metric:xic_fwhm")
})

test_that("ticQuantileRtFraction works properly.", {
    tmp <- ticQuantileRtFraction(chr)
    ## Data sorted by RT: (2.1,100), (2.5,250), (3.0,400), (3.4,300), (3.9,150),
    ##                    (5.1,80), (5.8,500), (6.3,1200), (6.9,600), (7.5,120)
    ## Cumulative TIC: 100, 350, 750, 1050, 1200, 1280, 1780, 2980, 3580, 3700
    ## Total TIC = 3700, RT duration = 5.4
    ## 0%: RT 2.1 → (2.1-2.1)/5.4 = 0
    ## 25%: cumsum >= 925 at RT 3.4 → (3.4-2.1)/5.4 = 0.2407
    ## 50%: cumsum >= 1850 at RT 6.3 → (6.3-2.1)/5.4 = 0.7778
    ## 75%: cumsum >= 2775 at RT 6.3 → (6.3-2.1)/5.4 = 0.7778
    ## 100%: cumsum >= 3700 at RT 7.5 → (7.5-2.1)/5.4 = 1.0
    expected <- c(0, (3.4-2.1)/5.4, (6.3-2.1)/5.4, (6.3-2.1)/5.4, 1.0)
    expect_equal(length(tmp), 5)
    expect_equal(names(tmp), c("0%", "25%", "50%", "75%", "100%"))
    expect_equal(as.numeric(tmp), expected, tolerance = 1e-4)

    ## test attributes
    expect_equal(attr(tmp, "ticQuantileRtFraction"), "MS:4000183")
})

test_that("areaUnderTic works properly.", {
    tmp <- areaUnderTic(chr)
    expect_equal(length(tmp), 1)
    expect_true(is.numeric(tmp))
    expect_equal(as.numeric(tmp), sum(c(100, 250, 400, 300, 150, 80, 500, 1200, 600, 120)))

    ## test attributes
    expect_equal(attr(tmp, "areaUnderTic"), "MS:4000155")
})

test_that("areaUnderTicRtQuantiles works properly.", {
    tmp <- areaUnderTicRtQuantiles(chr)
    expect_equal(length(tmp), 4)
    expect_equal(names(tmp), c("25%", "50%", "75%", "100%"))
    expect_true(all(is.numeric(tmp)))
    expect_true(sum(tmp) > 0)  # Total area should be positive

    ## test attributes
    expect_equal(attr(tmp, "areaUnderTicRtQuantiles"), "MS:4000156")
})

test_that("areaUnderTicRtQuantiles filters by msLevel for Chromatograms.", {
    cdata_mslevel <- data.frame(
        msLevel = c(1L, 2L),
        mz = c(100, 100),
        dataOrigin = c("mem1", "mem1")
    )
    pdata_mslevel <- list(
        data.frame(rtime = c(0, 1), intensity = c(1, 1)),
        data.frame(rtime = c(0, 1), intensity = c(2, 2))
    )
    chr_mslevel <- Chromatograms(ChromBackendMemory(), chromData = cdata_mslevel,
                                 peaksData = pdata_mslevel)

    res_ms1 <- areaUnderTicRtQuantiles(chr_mslevel, msLevel = 1L)
    expect_equal(as.numeric(res_ms1), rep(0.25, 4), tolerance = 1e-8)
    expect_equal(names(res_ms1), c("25%", "50%", "75%", "100%"))
    expect_equal(attr(res_ms1, "areaUnderTicRtQuantiles"), "MS:4000156")

    res_ms2 <- areaUnderTicRtQuantiles(chr_mslevel, msLevel = 2L)
    expect_equal(as.numeric(res_ms2), rep(0.5, 4), tolerance = 1e-8)
    expect_equal(names(res_ms2), c("25%", "50%", "75%", "100%"))

    res_empty <- areaUnderTicRtQuantiles(chr_mslevel, msLevel = 3L)
    expect_true(all(is.na(res_empty)))
    expect_equal(names(res_empty), c("25%", "50%", "75%", "100%"))
    expect_equal(attr(res_empty, "areaUnderTicRtQuantiles"), "MS:4000156")
})

## Chromatograms with MS1 and MS2 for areaUnderTicMs1/Ms2 tests
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
chr_ms <- Chromatograms(ChromBackendMemory(), chromData = cdata_ms, peaksData = pdata_ms)

test_that("areaUnderTicMs1 works properly.", {
    tmp <- areaUnderTicMs1(chr_ms)
    ## MS1 chromatograms: chr1 (100+200+300=600) + chr3 (400+500+600=1500) = 2100
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 2100)

    ## test attributes
    expect_equal(attr(tmp, "areaUnderTicMs1"), "MS:4000029")

    ## Test with no MS1 data
    cdata_ms2_only <- data.frame(msLevel = 2L, mz = 100.0, dataOrigin = "mem1")
    pdata_ms2_only <- list(data.frame(rtime = c(1, 2), intensity = c(10, 20)))
    chr_ms2_only <- Chromatograms(ChromBackendMemory(), chromData = cdata_ms2_only,
                                   peaksData = pdata_ms2_only)
    tmp2 <- areaUnderTicMs1(chr_ms2_only)
    expect_true(is.na(tmp2))
})

test_that("areaUnderTicMs2 works properly.", {
    tmp <- areaUnderTicMs2(chr_ms)
    ## MS2 chromatogram: chr2 (50+60+70=180)
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), 180)

    ## test attributes
    expect_equal(attr(tmp, "areaUnderTicMs2"), "MS:4000030")

    ## Test with no MS2 data
    cdata_ms1_only <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_ms1_only <- list(data.frame(rtime = c(1, 2), intensity = c(10, 20)))
    chr_ms1_only <- Chromatograms(ChromBackendMemory(), chromData = cdata_ms1_only,
                                   peaksData = pdata_ms1_only)
    tmp2 <- areaUnderTicMs2(chr_ms1_only)
    expect_true(is.na(tmp2))
})

