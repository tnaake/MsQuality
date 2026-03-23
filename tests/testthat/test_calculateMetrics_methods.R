## Test calculateMetrics S4 methods
fls <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
    MsDataHub::X20171016_POOL_POS_3_105.134.mzML())
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

test_that("areaUnderTic Chromatograms dispatch", {
    chr <- test_chr
    res <- areaUnderTic(chr)

    expect_equal(length(res), 2)
    expect_equal(res[1], sum(c(100, 250, 400, 300, 150)))
    expect_equal(res[2], sum(c(80, 500, 1200, 600, 120)))
    expect_equal(attr(res, "areaUnderTic"), "MS:4000155")
})

test_that("ticQuantileRtFraction Chromatograms dispatch", {
    chr <- test_chr
    res <- ticQuantileRtFraction(chr)

    expect_true(is.matrix(res))
    expect_equal(nrow(res), 2)
    expect_equal(ncol(res), 5)
    expect_equal(colnames(res), c("0%", "25%", "50%", "75%", "100%"))
    expect_true(all(res >= 0 & res <= 1, na.rm = TRUE))
    expect_equal(attr(res, "ticQuantileRtFraction"), "MS:4000183")
})

test_that("Chromatograms methods accept additional arguments and are consistent", {
    chr <- test_chr

    expect_no_error(areaUnderTic(chr))
    expect_no_error(ticQuantileRtFraction(chr, probs = c(0.25, 0.5, 0.75, 1)))

    custom <- ticQuantileRtFraction(chr, probs = c(0.5, 1.0))
    expect_true(is.matrix(custom))
    expect_equal(nrow(custom), 2)
    expect_equal(ncol(custom), 2)

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
    chrom_metrics <- c("chromatographyDuration", "maxIntensity", "intensityMean")

    result <- calculateMetrics(object = chr, metrics = chrom_metrics,
        filterEmptyObject = FALSE)

    expect_s3_class(result, "data.frame")
    expect_true(all(chrom_metrics[1] %in% colnames(result)))
    expect_true(is.numeric(result[, 1]))
})

test_that("calculateMetrics works with multiple Chromatograms metrics", {
    chr <- Chromatograms(spectra)
    chrom_metrics <- c("chromatographyDuration", "maxIntensity",
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

test_that("calculateMetrics returns correct column names with multiple metrics", {
    chr <- test_chr
    chrom_metrics <- c(
        "peakCount",
        "maxIntensity",
        "baselineIntensity",
        "signalToNoiseRatio"
    )

    result <- calculateMetrics(
        object = chr,
        metrics = chrom_metrics,
        filterEmptyObject = FALSE,
        na.rm = TRUE,
        probs = 0.05
    )

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_equal(colnames(result), chrom_metrics)
    expect_true(all(sapply(result, is.numeric)))
    expect_false(any(is.na(colnames(result))))
})

test_that("calculateMetrics works with metrics returning multiple values",
{    chr <- test_chr
    chrom_metrics <- c("xicFwhm", "gaussianSimilarity")
    expect_no_error({
        result <- calculateMetrics(
            object = chr,
            metrics = chrom_metrics,
            filterEmptyObject = FALSE
        )
    })

    result <- calculateMetrics(
        object = chr,
        metrics = chrom_metrics,
        filterEmptyObject = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_true("gaussianSimilarity.gaussian_similarity" %in% colnames(result))
    expect_true("gaussianSimilarity.gaussian_residuals" %in% colnames(result))
    expect_true("xicFwhm" %in% colnames(result))
})

test_that("calculateMetrics handles na.rm parameter without errors", {
    chr <- test_chr
    chrom_metrics <- c(
        "peakCount",
        "maxIntensity",
        "intensityMean",
        "intensitySd",
        "baselineIntensity",
        "signalToNoiseRatio"
    )

    expect_no_error({
        result <- calculateMetrics(
            object = chr,
            metrics = chrom_metrics,
            filterEmptyObject = FALSE,
            na.rm = TRUE,
            probs = 0.05
        )
    })

    result <- calculateMetrics(
        object = chr,
        metrics = chrom_metrics,
        filterEmptyObject = FALSE,
        na.rm = TRUE,
        probs = 0.05
    )

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), length(chrom_metrics))
    expect_equal(colnames(result), chrom_metrics)
})

test_that("calculateMetrics with signalToNoiseRatio does not conflict with method parameter", {

    chr <- test_chr
    chrom_metrics <- c("signalToNoiseRatio", "maxIntensity")

    ## Passing an unrelated 'method' parameter should not cause an error
    expect_no_error({
        result <- calculateMetrics(
            object = chr,
            metrics = chrom_metrics,
            filterEmptyObject = FALSE,
            method = "some_unrelated_value"
        )
    })

    result <- calculateMetrics(
        object = chr,
        metrics = chrom_metrics,
        filterEmptyObject = FALSE
    )
    expect_s3_class(result, "data.frame")
    expect_true("signalToNoiseRatio" %in% colnames(result))
})

################################################################################
######################### mzQC format for Chromatograms ########################
################################################################################

test_that("calculateMetricsFromChromatograms returns mzQC with PSI:MS metrics", {
    chr <- Chromatograms(spectra)

    suppressWarnings(
        res <- calculateMetricsFromChromatograms(
            chromatograms = chr,
            metrics = c("areaUnderTic", "ticQuantileRtFraction"),
            format = "mzQC"
        )
    )

    ## one MzQCmzQC object per sample
    expect_true(is.list(res))
    expect_equal(length(res), 2)
    expect_equal(class(res[[1]])[1], "MzQCmzQC")
    expect_equal(class(res[[2]])[1], "MzQCmzQC")

    ## each object should have exactly 2 quality metrics (both have MS: IDs)
    expect_equal(
        length(res[[1]]$runQualities[[1]]$qualityMetrics), 2)
    expect_equal(
        length(res[[2]]$runQualities[[1]]$qualityMetrics), 2)

    ## areaUnderTic
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession,
        "MS:4000155")
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name,
        "area under TIC")
    expect_true(
        is.numeric(res[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value))

    ## ticQuantileRtFraction
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[2]]$accession,
        "MS:4000183")
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[2]]$name,
        "TIC quantile RT fraction")
})

test_that("calculateMetricsFromChromatograms mzQC excludes custom metrics", {
    chr <- Chromatograms(spectra)

    ## maxIntensity has no PSI:MS attribute, so it should be
    ## excluded from the mzQC output
    suppressWarnings(
        res <- calculateMetricsFromChromatograms(
            chromatograms = chr,
            metrics = c("maxIntensity", "areaUnderTic"),
            format = "mzQC"
        )
    )

    ## only areaUnderTic should appear in the mzQC output
    expect_equal(
        length(res[[1]]$runQualities[[1]]$qualityMetrics), 1)
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession,
        "MS:4000155")
})

test_that("calculateMetrics generic dispatches mzQC for Chromatograms", {
    chr <- Chromatograms(spectra)

    suppressWarnings(
        res <- calculateMetrics(
            object = chr,
            metrics = c("areaUnderTic", "ticQuantileRtFraction"),
            format = "mzQC"
        )
    )

    expect_true(is.list(res))
    expect_equal(length(res), 2)
    expect_equal(class(res[[1]])[1], "MzQCmzQC")
    expect_equal(
        res[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession,
        "MS:4000155")
})

test_that("calculateMetricsFromChromatograms mzQC has URIs", {
    chr <- Chromatograms(spectra)

    suppressWarnings(
        res <- calculateMetricsFromChromatograms(
            chromatograms = chr,
            metrics = c("areaUnderTic"),
            format = "mzQC"
        )
    )

    expect_true(stringr::str_starts(
        res[[1]]$runQualities[[1]]$metadata$inputFiles[[1]]$location,
        "file://"))
})

test_that("calculateMetricsFromChromatograms mzQC structure is valid", {
    chr <- Chromatograms(spectra)

    suppressWarnings(
        res <- calculateMetricsFromChromatograms(
            chromatograms = chr,
            metrics = c("areaUnderTic"),
            format = "mzQC"
        )
    )

    ## check mzQC structure fields
    expect_true(!is.null(res[[1]]$version))
    expect_true(!is.null(res[[1]]$creationDate))
    expect_true(!is.null(res[[1]]$contactName))
    expect_true(!is.null(res[[1]]$description))
    expect_true(length(res[[1]]$runQualities) == 1)
    expect_true(length(res[[1]]$controlledVocabularies) == 1)

    ## metadata should contain analysisSoftware with MsQuality info
    software <- res[[1]]$runQualities[[1]]$metadata$analysisSoftware
    expect_true(length(software) >= 1)
})

test_that("calculateMetricsFromChromatograms data.frame format still works", {
    chr <- Chromatograms(spectra)

    res_df <- calculateMetricsFromChromatograms(
        chromatograms = chr,
        metrics = c("areaUnderTic", "maxIntensity"),
        format = "data.frame"
    )

    expect_s3_class(res_df, "data.frame")
    expect_equal(ncol(res_df), 2)
    expect_equal(colnames(res_df), c("areaUnderTic", "maxIntensity"))
})

test_that("mzQC from Chromatograms and Spectra agree on accession", {
    chr <- Chromatograms(spectra)

    suppressWarnings({
        res_chr <- calculateMetricsFromChromatograms(
            chromatograms = chr,
            metrics = c("areaUnderTic"),
            format = "mzQC"
        )
        res_sps <- calculateMetricsFromSpectra(
            spectra = spectra,
            metrics = c("areaUnderTic"),
            format = "mzQC",
            msLevel = 1
        )
    })

    ## both should produce the same accession and metric name
    expect_equal(
        res_chr[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession,
        res_sps[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession)
    expect_equal(
        res_chr[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name,
        res_sps[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name)

    ## both values should be numeric
    expect_true(
        is.numeric(res_chr[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value))
    expect_true(
        is.numeric(res_sps[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value))
})
