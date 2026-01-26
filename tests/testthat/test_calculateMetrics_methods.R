## Test calculateMetrics S4 methods
library("Spectra")
library("MsExperiment")
library("Chromatograms")
library("S4Vectors")
library("MsQuality")

fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## Shared Chromatograms fixture
test_chr <- Chromatograms(
    chromData = data.frame(
        msLevel = c(1L, 1L),
        mz = c(112.2, 123.3),
        dataOrigin = c("mem1", "mem1")
    ),
    peaksData = list(
        data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
            intensity = c(100, 250, 400, 300, 150)),
        data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
            intensity = c(80, 500, 1200, 600, 120))
    ),
    backend = ChromBackendMemory()
)

## Define quality metrics to test
metrics <- c("chromatographyDuration", "ticQuantileRtFraction", "numberSpectra")

test_that("areaUnderTic Chromatograms method dispatch", {
    expect_true(hasMethod("areaUnderTic", "Chromatograms"))

    chr <- test_chr
    res <- areaUnderTic(chr)

    expect_equal(as.numeric(res), 3700)
    expect_equal(attr(res, "areaUnderTic"), "MS:4000155")
})

test_that("ticQuantileRtFraction Chromatograms method dispatch", {
    expect_true(hasMethod("ticQuantileRtFraction", "Chromatograms"))

    chr <- test_chr
    res <- ticQuantileRtFraction(chr)

    expect_length(res, 5)
    expect_equal(names(res), c("0%", "25%", "50%", "75%", "100%"))
    expect_true(all(res >= 0 & res <= 1, na.rm = TRUE))
    expect_equal(attr(res, "ticQuantileRtFraction"), "MS:4000183")
})

test_that("Chromatograms methods accept additional arguments and are consistent", {
    chr <- test_chr

    expect_no_error(areaUnderTic(chr))
    expect_no_error(ticQuantileRtFraction(chr, probs = c(0.25, 0.5, 0.75, 1)))

    custom <- ticQuantileRtFraction(chr, probs = c(0.5, 1.0))
    expect_length(custom, 2)

    expect_equal(areaUnderTic(chr), areaUnderTic(chr))
    expect_equal(ticQuantileRtFraction(chr), ticQuantileRtFraction(chr))
})

test_that("calculateMetrics works with Spectra objects", {
    result <- calculateMetrics(object = spectra, metrics = metrics,
        filterEmptyObject = FALSE, msLevel = 1)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_true(all(metrics[1] %in% colnames(result)))
})

test_that("calculateMetrics works with multiple Spectra metrics", {
    spectra_metrics <- c("chromatographyDuration", "numberSpectra", "rtAcquisitionRange")

    result <- calculateMetrics(object = spectra, metrics = spectra_metrics,
        filterEmptyObject = FALSE, msLevel = 1)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    for (metric in spectra_metrics) {
        expect_true(any(grepl(metric, colnames(result), fixed = TRUE)),
                   info = paste("Metric", metric, "not found in result"))
    }
})

test_that("calculateMetrics works with MsExperiment objects", {
    msexp <- MsExperiment()
    sd <- DataFrame(sample_id = c("QC1", "QC2"),
        sample_name = c("QC Pool", "QC Pool"))
    sampleData(msexp) <- sd
    experimentFiles(msexp) <- MsExperimentFiles(mzML_files = fls)
    msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
        sampleIndex = c(1, 2), withIndex = c(1, 2))
    spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())

    result <- calculateMetrics(object = msexp, metrics = metrics,
        filterEmptyObject = FALSE, msLevel = 1)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_true(all(metrics[1] %in% colnames(result)))
})

test_that("calculateMetrics works with Chromatograms objects", {
    chr <- Chromatograms(spectra)
    chrom_metrics <- c("chromatogramDuration", "maxIntensity", "intensityMean")

    result <- calculateMetrics(object = chr, metrics = chrom_metrics,
        filterEmptyObject = FALSE)

    expect_s3_class(result, "data.frame")
    expect_true(all(chrom_metrics[1] %in% colnames(result)))
    expect_true(is.numeric(result[, 1]))
})

test_that("calculateMetrics works with multiple Chromatograms metrics", {
    chr <- Chromatograms(spectra)
    chrom_metrics <- c("chromatogramDuration", "maxIntensity",
                       "intensityMean", "intensitySd", "peakCount")

    result <- calculateMetrics(object = chr, metrics = chrom_metrics,
        filterEmptyObject = FALSE)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    for (metric in chrom_metrics) {
        expect_true(any(grepl(metric, colnames(result), fixed = TRUE)),
                   info = paste("Metric", metric, "not found in result"))
    }
    for (col in chrom_metrics) {
        col_idx <- which(grepl(col, colnames(result), fixed = TRUE))[1]
        expect_true(is.numeric(result[, col_idx]),
                   info = paste("Column", col, "is not numeric"))
    }
})

test_that("calculateMetrics dispatches correctly based on class", {
    result_spectra <- calculateMetrics(spectra, metrics = metrics, msLevel = 1)
    expect_s3_class(result_spectra, "data.frame")

    msexp <- MsExperiment()
    spectra(msexp) <- spectra
    result_msexp <- calculateMetrics(msexp, metrics = metrics, msLevel = 1)
    expect_s3_class(result_msexp, "data.frame")

    expect_equal(result_spectra, result_msexp)
})

test_that("calculateMetrics handles metrics with different parameter requirements", {
    spectra_metrics <- c("numberSpectra", "rtAcquisitionRange", "chromatographyDuration")

    result_ms1 <- calculateMetrics(object = spectra, metrics = spectra_metrics,
        filterEmptyObject = FALSE, msLevel = 1)
    result_ms2 <- calculateMetrics(object = spectra, metrics = spectra_metrics,
        filterEmptyObject = FALSE, msLevel = 2)

    expect_s3_class(result_ms1, "data.frame")
    expect_s3_class(result_ms2, "data.frame")
    expect_equal(nrow(result_ms1), 2)
    expect_equal(nrow(result_ms2), 2)
})

test_that("calculateMetrics passes additional parameters to metric functions", {
    spectra_metrics <- c("msSignal10xChange")

    result_jump <- calculateMetrics(object = spectra, metrics = spectra_metrics,
        filterEmptyObject = FALSE, msLevel = 1, change = "jump")
    result_fall <- calculateMetrics(object = spectra, metrics = spectra_metrics,
        filterEmptyObject = FALSE, msLevel = 1, change = "fall")

    expect_s3_class(result_jump, "data.frame")
    expect_s3_class(result_fall, "data.frame")
    expect_equal(nrow(result_jump), 2)
    expect_equal(nrow(result_fall), 2)
    expect_true(all(is.numeric(as.matrix(result_jump))))
    expect_true(all(is.numeric(as.matrix(result_fall))))
})

test_that("calculateMetrics validates inputs correctly", {
    expect_error(
        calculateMetrics(spectra, metrics = "invalid_metric"),
        "should be one of"
    )

    expect_error(
        calculateMetrics(spectra, metrics = metrics, filterEmptyObject = "yes"),
        "'filterEmptyObject' has to be either TRUE or FALSE"
    )

    expect_error(
        calculateMetrics(Chromatograms(spectra), metrics = "invalid_metric"),
        "should be one of"
    )
})
