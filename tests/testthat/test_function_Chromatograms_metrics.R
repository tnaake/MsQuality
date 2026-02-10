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
    ## Returns total count of data points across all chromatograms
    tmp <- peakCount(chr)
    expect_equal(length(tmp), 1)        # Single aggregated value
    expect_equal(as.numeric(tmp), 10)   # 5 + 0 + 5 = 10 total points

    ## test attributes
    expect_equal(attr(tmp, "peakCount"), "custom_metric:peak_count")
})

test_that("peakCount works with na.rm parameter.", {
    ## Test with all non-NA data (standard chr object)
    ## na.rm = TRUE and na.rm = FALSE should give same result when no NAs
    tmp_standard_rm <- peakCount(chr, na.rm = TRUE)
    expect_equal(as.numeric(tmp_standard_rm), 10)  # All 10 points are non-NA

    tmp_standard_no_rm <- peakCount(chr, na.rm = FALSE)
    expect_equal(as.numeric(tmp_standard_no_rm), 10)  # Same result when no NAs
    expect_equal(attr(tmp_standard_no_rm, "peakCount"), "custom_metric:peak_count")
})

test_that("rtIqr works properly for Chromatograms.", {
    tmp <- rtIqr(chr)
    ## Aggregated IQR across all chromatograms
    expect_equal(length(tmp), 1)
    expect_equal(as.numeric(tmp), IQR(c(2.1, 2.5, 3.0, 3.4, 3.9, 5.1, 5.8, 6.3, 6.9, 7.5)))
    ## test attributes
    expect_equal(attr(tmp, "rtIqr"), "custom_metric:rt_iqr")
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

test_that("xicFwhm works properly.", {
    ## xicFwhm is for single chromatograms (EIC/XIC) - use chr[1]
    tmp <- xicFwhm(chr[1])
    expect_equal(length(tmp), 1)
    expect_true(is.numeric(tmp))
    expect_true(tmp > 0)  # FWHM should be positive

    ## Test with empty chromatogram returns NA
    tmp_empty <- xicFwhm(chr[2])
    expect_true(is.na(tmp_empty))

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

## peakBoundary tests
test_that("peakBoundary returns correct RT boundaries for clean symmetric peak.", {
    ## Simple chromatogram with clear symmetric peak dropping to zero
    cdata_pb <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pb <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb, peaksData = pdata_pb)

    ## Uses MsCoreUtils::valleys() to find local minima
    tmp <- peakBoundary(chr_pb)
    expect_equal(length(tmp), 2)
    expect_equal(names(tmp), c("peakBoundary_left", "peakBoundary_right"))
    ## Peak apex at rt=4, valleys() finds valleys at the zeros (rt 1 and 7)
    expect_equal(unname(tmp["peakBoundary_left"]), 1)
    expect_equal(unname(tmp["peakBoundary_right"]), 7)
    expect_equal(attr(tmp, "peakBoundary"), "custom_metric:peak_boundary")
})

test_that("peakBoundary finds first valley in overlapping peaks.", {
    ## Two overlapping peaks with a valley between them
    cdata_pb <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pb <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
        intensity = c(10, 50, 100, 60, 30, 60, 100, 50, 10)  # valley at rt=5
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb, peaksData = pdata_pb)

    tmp <- peakBoundary(chr_pb)  # finds highest peak (rt=3 or rt=7)
    expect_equal(length(tmp), 2)
    ## valleys() should stop at the valley between peaks
    expect_true(tmp["peakBoundary_left"] >= 1)
    expect_true(tmp["peakBoundary_right"] <= 9)
})

test_that("peakBoundary handles peak with elevated baseline.", {
    ## Peak on elevated baseline - should find where signal levels off
    cdata_pb <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pb <- list(data.frame(
        rtime = 1:15,
        intensity = c(100, 100, 100, 150, 300, 500, 300, 150, 100, 100, 100, 100, 100, 100, 100)
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb, peaksData = pdata_pb)

    tmp <- peakBoundary(chr_pb)
    expect_equal(length(tmp), 2)
    ## Peak apex at rt=6, should find where it becomes flat (around rt 3 and rt 9)
    expect_true(tmp["peakBoundary_left"] >= 1 && tmp["peakBoundary_left"] <= 4)
    expect_true(tmp["peakBoundary_right"] >= 8 && tmp["peakBoundary_right"] <= 11)
})

test_that("peakBoundary handles tailing peak.", {
    ## Asymmetric peak: sharp rise, slow decay (common in chromatography)
    cdata_pb <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pb <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        intensity = c(5, 20, 100, 80, 50, 35, 25, 20, 18, 15)
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb, peaksData = pdata_pb)

    tmp <- peakBoundary(chr_pb)
    expect_equal(length(tmp), 2)
    ## Peak apex at rt=3, left boundary should be at first point
    expect_equal(unname(tmp["peakBoundary_left"]), 1)
    ## Right boundary extends to end since no valley (monotonic decrease)
    expect_equal(unname(tmp["peakBoundary_right"]), 10)
})

test_that("peakBoundary handles real-world peak with elevated baseline and noisy tail.", {
    ## Real-world example: peak with high baseline and long noisy tail
    ## Peak apex around index 12, then noisy plateau that never drops to baseline
    cdata_pb <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pb <- list(data.frame(
        rtime = seq(1, 148, by = 1),
        intensity = c(
            3011.207, 2974.531, 2762.167, 2891.680, 3585.870, 3495.168, 4962.872,
            14508.343, 45776.547, 99684.969, 144267.156, 145303.016, 117992.516,
            81776.422, 56023.023, 43023.773, 35532.582, 32609.395, 27691.934,
            25894.314, 25407.229, 27242.326, 24835.045, 23004.254, 20529.389,
            22985.701, 19728.000, 21125.123, 17173.430, 16828.990, 17671.371,
            17238.762, 14417.884, 16073.300, 15524.321, 15545.252, 16164.231,
            14967.193, 14825.901, 14622.387, 14367.496, 14099.743, 13814.941,
            13847.837, 13400.929, 14825.461, 14537.347, 14367.999, 15478.081,
            13988.006, 15239.449, 13623.584, 15051.759, 13910.513, 12910.956,
            13675.151, 13019.853, 12531.170, 13472.155, 11975.559, 12469.830,
            14894.729, 13566.914, 12935.803, 15064.299, 13708.584, 15143.816,
            12255.616, 13121.985, 13695.070, 12550.412, 13944.659, 13570.368,
            13099.557, 13503.188, 10923.682, 12251.983, 12864.445, 12905.133,
            13073.855, 12577.082, 13440.513, 15015.454, 14148.987, 14335.039,
            13118.187, 14612.891, 12529.094, 13824.206, 14526.008, 15842.500,
            15751.920, 15425.752, 15308.919, 14428.304, 14696.171, 14999.367,
            13332.374, 14872.715, 12422.044, 13559.954, 12324.938, 15466.142,
            15311.044, 13859.441, 14940.694, 15187.022, 14833.513, 15504.637,
            14668.642, 13930.903, 14814.932, 14515.681, 14756.643, 15067.479,
            14381.385, 16005.971, 14378.000, 16044.934, 16053.924, 15861.024,
            15804.795, 15378.197, 16378.093, 15469.408, 16802.516, 15320.836,
            14691.501, 14563.201, 15646.632, 14151.453, 16516.934, 15049.130,
            14293.031, 14485.947, 14089.717, 15332.651, 17004.359, 13967.167,
            14529.343, 13746.593, 15503.060, 13467.515, 14113.591, 14057.125,
            12294.305, 13497.132, 12305.407
        )
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb, peaksData = pdata_pb)

    tmp <- peakBoundary(chr_pb)
    expect_equal(length(tmp), 2)
    expect_equal(names(tmp), c("peakBoundary_left", "peakBoundary_right"))

    ## Peak apex is at index 12 (rt=12), max intensity = 145303
    ## MsCoreUtils::valleys() finds boundaries at indices 6 and 21
    ## This is where the noisy baseline starts/ends (local minima)
    expect_equal(unname(tmp["peakBoundary_left"]), 6)
    expect_equal(unname(tmp["peakBoundary_right"]), 21)

    ## The peak width should be reasonable (not the entire chromatogram)
    peak_width <- unname(tmp["peakBoundary_right"] - tmp["peakBoundary_left"])
    expect_equal(peak_width, 15)  # 21 - 6 = 15
})

test_that("peakBoundary returns NA for empty chromatograms.", {
    cdata_empty <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_empty <- list(data.frame(rtime = numeric(), intensity = numeric()))
    chr_empty <- Chromatograms(ChromBackendMemory(), chromData = cdata_empty, peaksData = pdata_empty)

    tmp <- peakBoundary(chr_empty)
    expect_equal(length(tmp), 2)
    expect_true(all(is.na(tmp)))
})

test_that("peakBoundary handles max intensity at first position.", {
    cdata_first <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_first <- list(data.frame(
        rtime = seq(1, 18, length.out = 18),
        intensity = c(4174.2911, 2450.8100, 530.3057, 963.4568, 1058.0544, 
                      2244.1909, 3708.5597, 4036.2473, 2960.1146, 3246.5288, 
                      2472.2660, 1355.1185, 2020.3214, 1210.5110, 718.1166, 
                      900.8450, 799.6667, 1926.9595)
    ))
    chr_first <- Chromatograms(ChromBackendMemory(), chromData = cdata_first, 
                               peaksData = pdata_first)

    tmp <- peakBoundary(chr_first)
    expect_equal(length(tmp), 2)
    expect_false(any(is.na(tmp)))
    ## Left boundary should be at first position (index 1)
    expect_equal(unname(tmp["peakBoundary_left"]), 1)
})

test_that("peakBoundary handles max intensity at last position.", {
    cdata_last <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_last <- list(data.frame(
        rtime = seq(1, 10, length.out = 10),
        intensity = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
    ))
    chr_last <- Chromatograms(ChromBackendMemory(), chromData = cdata_last, 
                              peaksData = pdata_last)

    tmp <- peakBoundary(chr_last)
    expect_equal(length(tmp), 2)
    expect_false(any(is.na(tmp)))
    ## Right boundary should be at last position
    expect_equal(unname(tmp["peakBoundary_right"]), 10)
})

## peakWidth tests
test_that("peakWidth returns correct width for clean peak.", {
    cdata_pw <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pw <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_pw <- Chromatograms(ChromBackendMemory(), chromData = cdata_pw, peaksData = pdata_pw)

    tmp <- peakWidth(chr_pw)
    expect_equal(length(tmp), 1)
    ## valleys() finds valleys at rt 1 and 7, width = 6
    expect_equal(as.numeric(tmp), 6)
    expect_equal(attr(tmp, "peakWidth"), "custom_metric:peak_width")
})

test_that("peakWidth returns NA for empty chromatograms.", {
    cdata_empty <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_empty <- list(data.frame(rtime = numeric(), intensity = numeric()))
    chr_empty <- Chromatograms(ChromBackendMemory(), chromData = cdata_empty, peaksData = pdata_empty)

    tmp <- peakWidth(chr_empty)
    expect_equal(length(tmp), 1)
    expect_true(is.na(tmp))
})

test_that("peakWidth accepts pre-computed peakBoundary.", {
    cdata_pw <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_pw <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_pw <- Chromatograms(ChromBackendMemory(), chromData = cdata_pw, peaksData = pdata_pw)

    ## Calculate peakBoundary once
    pb <- peakBoundary(chr_pw)

    ## Pass to peakWidth
    tmp <- peakWidth(chr_pw, peakBoundary = pb)
    expect_equal(as.numeric(tmp), 6)
})

## peakBeta tests
test_that("peakBeta returns beta values for a clean peak.", {
    ## Create a bell-shaped peak with enough points
    cdata_beta <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_beta <- list(data.frame(
        rtime = seq(1, 20, by = 0.5),
        intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
                      700, 650, 500, 350, 200, 100, 50, 25, 15,
                      10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1)
    ))
    chr_beta <- Chromatograms(ChromBackendMemory(), chromData = cdata_beta, peaksData = pdata_beta)

    tmp <- peakBeta(chr_beta)
    expect_equal(length(tmp), 2)
    expect_equal(names(tmp), c("beta_cor", "beta_snr"))
    expect_true(is.numeric(tmp["beta_cor"]))
    expect_true(is.numeric(tmp["beta_snr"]))
    expect_equal(attr(tmp, "peakBeta"), "custom_metric:peak_beta")
})

test_that("peakBeta accepts pre-computed peakBoundary.", {
    cdata_beta <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_beta <- list(data.frame(
        rtime = seq(1, 20, by = 0.5),
        intensity = c(10, 15, 25, 50, 100, 200, 350, 500, 650,
                      700, 650, 500, 350, 200, 100, 50, 25, 15,
                      10, 8, 6, 5, 4, 3, 2, 2, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1)
    ))
    chr_beta <- Chromatograms(ChromBackendMemory(), chromData = cdata_beta, peaksData = pdata_beta)

    pb <- peakBoundary(chr_beta)
    tmp <- peakBeta(chr_beta, peakBoundary = pb)
    expect_equal(length(tmp), 2)
    expect_true(!is.na(tmp[["beta_cor"]]))
})

test_that("peakBeta returns NA for empty chromatograms.", {
    cdata_empty <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_empty <- list(data.frame(rtime = numeric(), intensity = numeric()))
    chr_empty <- Chromatograms(ChromBackendMemory(), chromData = cdata_empty, peaksData = pdata_empty)

    tmp <- peakBeta(chr_empty)
    expect_equal(length(tmp), 2)
    expect_true(all(is.na(tmp)))
})

test_that("peakBeta returns NA when peak has < 5 points.", {
    ## Very short peak
    cdata_short <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_short <- list(data.frame(
        rtime = c(1, 2, 3),
        intensity = c(10, 100, 10)
    ))
    chr_short <- Chromatograms(ChromBackendMemory(), chromData = cdata_short, peaksData = pdata_short)

    tmp <- peakBeta(chr_short)
    expect_equal(length(tmp), 2)
    expect_true(all(is.na(tmp)))
})

test_that("xicFwhm accepts pre-computed peakBoundary.", {
    cdata_fwhm <- data.frame(msLevel = 1L, mz = 100.0, dataOrigin = "mem1")
    pdata_fwhm <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_fwhm <- Chromatograms(ChromBackendMemory(), chromData = cdata_fwhm, peaksData = pdata_fwhm)

    pb <- peakBoundary(chr_fwhm)
    tmp <- xicFwhm(chr_fwhm, peakBoundary = pb)
    expect_equal(length(tmp), 1)
    expect_true(!is.na(tmp))
    expect_true(tmp > 0)
})

test_that("metrics handle data with NA values correctly.", {
    ## Create chromatogram with NA values in intensity
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

    ## With na.rm = TRUE, should compute without error and return valid numeric
    tmp_max <- maxIntensity(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_max), 500)

    tmp_mean <- intensityMean(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_mean), mean(c(100, 400, 150, 80, 500), na.rm = TRUE))

    tmp_sd <- intensitySd(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_sd), sd(c(100, 400, 150, 80, 500), na.rm = TRUE))

    tmp_range <- intensityRange(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_range), c(80, 500))

    tmp_baseline <- baselineIntensity(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_baseline),
                 as.numeric(quantile(c(100, 400, 150, 80, 500), probs = 0.05, na.rm = TRUE)))

    ## peakCount with na.rm should count only non-NA values
    tmp_count <- peakCount(chr_na, na.rm = TRUE)
    expect_equal(as.numeric(tmp_count), 5)  # 5 non-NA values

    tmp_count_all <- peakCount(chr_na, na.rm = FALSE)
    expect_equal(as.numeric(tmp_count_all), 8)  # 8 total values including NAs
})

