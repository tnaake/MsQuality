################################################################################
###################### format = 'data.frame' (default) #########################
################################################################################

fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## build the results
## define the quality metrics to be calculated
metrics <- c("chromatographyDuration", "ticQuartersRtFraction", 
    "ticQuartileToQuartileLogRatio", "numberSpectra", "areaUnderTic", 
    "msSignal10xChange")

## additional parameters passed to the quality metrics functions
## (msLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
suppressWarnings(
    metrics_spectra <- calculateMetricsFromSpectra(spectra = spectra,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_spectra_filtered <- calculateMetricsFromSpectra(spectra = spectra,
        metrics = metrics, filterEmptySpectra = TRUE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))

## calculate the metrics from Spectra
dO <- unique(spectra$dataOrigin)
spectra_1 <- spectra[spectra$dataOrigin == dO[1], ]
spectra_2 <- spectra[spectra$dataOrigin == dO[2], ]
suppressWarnings(metrics_spectra_1 <- calculateMetricsFromOneSampleSpectra(
    spectra = spectra_1, metrics = metrics, filterEmptySpectra = FALSE, 
    msLevel = 1, relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_spectra_1_filtered <- calculateMetricsFromOneSampleSpectra(
        spectra = spectra_1, metrics = metrics, filterEmptySpectra = TRUE, 
        msLevel = 1, relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(metrics_spectra_2 <- calculateMetricsFromOneSampleSpectra(
    spectra = spectra_2, metrics = metrics, filterEmptySpectra = FALSE, 
    msLevel = 1, relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_spectra_2_filtered <- calculateMetricsFromOneSampleSpectra(
        spectra = spectra_2, metrics = metrics, filterEmptySpectra = TRUE, 
        msLevel = 1, relativeTo = "Q1", mode = "TIC", change = "jump"))

## START unit test calculateMetricsFromOneSampleSpectra ## 
colnames_metrics <- c("chromatographyDuration", "ticQuartersRtFraction.0%",                  
    "ticQuartersRtFraction.25%", "ticQuartersRtFraction.50%",
    "ticQuartersRtFraction.75%", "ticQuartersRtFraction.100%",
    "ticQuartileToQuartileLogRatio.Q2/Q1",
    "ticQuartileToQuartileLogRatio.Q3/Q1",
    "ticQuartileToQuartileLogRatio.Q4/Q1",
    "numberSpectra", "areaUnderTic", "msSignal10xChange")
metrics_spectra_1_vals <- c(2.594770e+02, 0, 2.505386e-01,
    4.999981e-01, 7.505367e-01, 1.000000e+00, -1.603689e-01, -7.040791e-01,
    -5.018226e-01, 9.310000e02, 6.509697e+08, 0)
metrics_spectra_2_vals <- c(2.594770e+02, 0, 2.505386e-01,
    4.999981e-01, 7.505367e-01, 1.000000e+00, 5.052683e-02, -5.673914e-01,
    -4.700149e-01, 9.310000e02, 6.229579e+08, 0)

## create small test spectra to test filterEmptySpectra
spd <- DataFrame(
    msLevel = c(2L, 2L, 2L),
    polarity = c(1L, 1L, 1L),
    name = c("comp_1", "comp_2", "comp_3"))
## Assign m/z and intensity values
spd$mz <- list(
        c(109.2),
        c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
        c())
spd$intensity <- list(
    c(0),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c())
spd$rtime <- c(9.44, 9.44, 15.84)
sps_empty <- Spectra(spd)
dataOrigin(sps_empty) <- c("manual", "manual", "manual")
spd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
    c())
spd$intensity <- list(
    c(0, 0, 0, 0, 0),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c())
spd$rtime <- c(9.44, 9.44, 15.84)
sps_multiple_empty <- Spectra(spd)
dataOrigin(sps_multiple_empty) <- c("manual", "manual", "manual")
spd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
    c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
        111.0551, 123.0429, 138.0662, 195.0876))
spd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
spd$rtime <- c(9.44, 9.44, 15.84)
sps_not_empty <- Spectra(spd)
dataOrigin(sps_not_empty) <- c("manual", "manual", "manual")

test_that("calculateMetricsFromOneSampleSpectra", {
    ## spectra
    expect_error(
        calculateMetricsFromOneSampleSpectra(spectra, metrics = metrics),
        "'spectra' should only contain data from one origin")
    
    ## spectra_1
    expect_true(is.numeric(metrics_spectra_1))
    expect_equal(length(metrics_spectra_1), 12)
    expect_equal(names(metrics_spectra_1), colnames_metrics)
    expect_true(is.numeric(metrics_spectra_1))
    expect_equal(as.numeric(metrics_spectra_1), metrics_spectra_1_vals, 
        tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_1_filtered), metrics_spectra_1_vals, 
        tolerance = 1e-06)
    expect_error(calculateMetricsFromOneSampleSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromOneSampleSpectra("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_1, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(
        calculateMetricsFromOneSampleSpectra(spectra_1, metrics = "foo"),
        "should be one of")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_1, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
    
    ## spectra_2
    expect_true(is.numeric(metrics_spectra_2))
    expect_true(is.numeric(metrics_spectra_2_filtered))
    expect_equal(length(metrics_spectra_2), 12)
    expect_equal(length(metrics_spectra_2_filtered), 12)
    expect_equal(names(metrics_spectra_2), colnames_metrics)
    expect_equal(names(metrics_spectra_2_filtered), colnames_metrics)
    expect_true(is.numeric(metrics_spectra_2))
    expect_true(is.numeric(metrics_spectra_2_filtered))
    expect_equal(as.numeric(metrics_spectra_2), metrics_spectra_2_vals, 
        tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_2_filtered), metrics_spectra_2_vals, 
        tolerance = 1e-06)
    expect_error(calculateMetricsFromOneSampleSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_2, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_2,
        metrics = "msSignal10xChange", change = c("jump", "fall")),
        "'change' has to be of length 1")
    
    ## test filterEmptySpectra
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_multiple_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_multiple_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_not_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromOneSampleSpectra(
        spectra = sps_not_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 3)
    
    ## test attributes
    expect_equal(attr(metrics_spectra_1, "names"), colnames_metrics)
    expect_equal(attr(metrics_spectra_1_filtered, "names"), colnames_metrics)
    expect_equal(attr(metrics_spectra_1, "chromatographyDuration"), "MS:4000053")
    expect_equal(attr(metrics_spectra_1_filtered, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra_1, "ticQuartersRtFraction"), "MS:4000054")
    expect_equal(attr(metrics_spectra_1_filtered, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_spectra_1, "numberSpectra"), "MS:4000059")
    expect_equal(attr(metrics_spectra_1_filtered, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_spectra_1, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_spectra_1_filtered, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_spectra_1, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_1_filtered, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_1, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_1_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_1, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_1_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_1, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_1_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_1, "change"), "jump")
    expect_equal(attr(metrics_spectra_1_filtered, "change"), "jump")
    expect_equal(attr(metrics_spectra_2, "names"), colnames_metrics)
    expect_equal(attr(metrics_spectra_2_filtered, "names"), colnames_metrics)
    expect_equal(attr(metrics_spectra_2, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra_2_filtered, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra_2, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_spectra_2_filtered, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_spectra_2, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_spectra_2_filtered, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_spectra_2, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_spectra_2_filtered, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_spectra_2, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_2_filtered, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_2, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_2_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_2, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_2_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_2, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_2_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_2, "change"), "jump")
    expect_equal(attr(metrics_spectra_2_filtered, "change"), "jump")
})
## END unit test calculateMetricsFromOneSampleSpectra

## START unit test calculateMetricsFromSpectra ##
test_that("calculateMetricsFromSpectra", {
    expect_equal(dim(metrics_spectra), c(2, 12))
    expect_equal(dim(metrics_spectra_filtered), c(2, 12))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_spectra), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
          "20171016_POOL_POS_3_105-134.mzML"))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_spectra_filtered), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    expect_equal(colnames(metrics_spectra), colnames_metrics)
    expect_equal(colnames(metrics_spectra_filtered), colnames_metrics)
    expect_equal(as.numeric(metrics_spectra[1, ]), 
        metrics_spectra_1_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_filtered[1, ]), 
        metrics_spectra_1_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra[2, ]), 
        metrics_spectra_2_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_filtered[2, ]), 
        metrics_spectra_2_vals, tolerance = 1e-06)
    
    expect_error(calculateMetricsFromSpectra("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromSpectra(spectra, metrics = "foo"),
        "should be one of ")
    expect_error(calculateMetricsFromSpectra(spectra, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(calculateMetricsFromSpectra(spectra, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
    
    ## test filterEmptySpectra
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_multiple_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_multiple_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_not_empty, metrics = "numberSpectra", 
        filterEmptySpectra = FALSE, msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetricsFromSpectra(
        spectra = sps_not_empty, metrics = "numberSpectra", 
        filterEmptySpectra = TRUE, msLevel = 2L)), 3)
    
    ## test attributes
    expect_equal(attr(metrics_spectra, "names"), NULL)
    expect_equal(attr(metrics_spectra_filtered, "names"), NULL)
    expect_equal(attr(metrics_spectra, "chromatographyDuration"), "MS:4000053")
    expect_equal(attr(metrics_spectra_filtered, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra, "ticQuartersRtFraction"), "MS:4000054")
    expect_equal(attr(metrics_spectra_filtered, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_spectra, "numberSpectra"), "MS:4000059")
    expect_equal(attr(metrics_spectra_filtered, "numberSpectra"), "MS:4000059")
    expect_equal(attr(metrics_spectra, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_spectra_filtered, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_spectra, "msSignal10xChange"), "MS:4000097")
    expect_equal(attr(metrics_spectra_filtered, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_spectra, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_spectra, "change"), "jump")
    expect_equal(attr(metrics_spectra_filtered, "change"), "jump")
})
## END unit test calculateMetricsFromSpectra ##


## START unit test calculateMetricsFromMsExperiment ##
msexp <- MsExperiment()
sd <- DataFrame(sample_id = c("QC1", "QC2"),
     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
sampleData(msexp) <- sd
 
## define file names containing spectra data for the samples and
## add them, along with other arbitrary files to the experiment
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
experimentFiles(msexp) <- MsExperimentFiles(
    mzML_files = fls,
    annotations = "internal_standards.txt")
## link samples to data files: first sample to first file in "mzML_files",
## second sample to second file in "mzML_files"
msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
    sampleIndex = c(1, 2), withIndex = c(1, 2))
msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
    sampleIndex = c(1, 2), withIndex = c(1, 1))

## import the data and add it to the mse object
spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())
 
## additional parameters passed to the quality metrics functions
## (msLevel is an argument of areaUnderTic and msSignal10xChange,
## relativeTo is an argument of msSignal10xChange) passed to ...
metrics_msexp <- calculateMetricsFromMsExperiment(msexp = msexp, 
    metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
    relativeTo = "Q1", mode = "TIC", change = "jump")
metrics_msexp_filtered <- calculateMetricsFromMsExperiment(msexp = msexp, 
    metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
    relativeTo = "Q1", mode = "TIC", change = "jump")

test_that("calculateMetricsFromMsExperiment", {
    expect_equal(dim(metrics_msexp), c(2, 12))
    expect_equal(dim(metrics_msexp_filtered), c(2, 12))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_msexp), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_msexp_filtered), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    expect_equal(colnames(metrics_msexp), colnames_metrics)
    expect_equal(colnames(metrics_msexp_filtered), colnames_metrics)
    expect_equal(as.numeric(metrics_msexp[1, ]), 
        metrics_spectra_1_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_filtered[1, ]), 
        metrics_spectra_1_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp[2, ]), 
        metrics_spectra_2_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_filtered[2, ]), 
        metrics_spectra_2_vals, tolerance = 1e-06)
    
    expect_error(calculateMetricsFromMsExperiment("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromMsExperiment(msexp, metrics = "foo"),
        "should be one of ")
    expect_error(calculateMetricsFromMsExperiment(msexp, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(calculateMetricsFromMsExperiment(msexp, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
    
    ## test attributes
    expect_equal(attr(metrics_msexp, "names"), NULL)
    expect_equal(attr(metrics_msexp_filtered, "names"), NULL)
    expect_equal(attr(metrics_msexp, "chromatographyDuration"), "MS:4000053")
    expect_equal(attr(metrics_msexp_filtered, "chromatographyDuration"), 
                 "MS:4000053")
    expect_equal(attr(metrics_msexp, "ticQuartersRtFraction"), "MS:4000054")
    expect_equal(attr(metrics_msexp_filtered, "ticQuartersRtFraction"), 
                 "MS:4000054")
    expect_equal(attr(metrics_msexp, "numberSpectra"), "MS:4000059")
    expect_equal(attr(metrics_msexp_filtered, "numberSpectra"), "MS:4000059")
    expect_equal(attr(metrics_msexp, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_msexp_filtered, "areaUnderTic"), "MS:4000155")
    expect_equal(attr(metrics_msexp, "msSignal10xChange"), "MS:4000097")
    expect_equal(attr(metrics_msexp_filtered, "msSignal10xChange"), 
                 "MS:4000097")
    expect_equal(attr(metrics_msexp, "msLevel"), 1)
    expect_equal(attr(metrics_msexp_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_msexp, "relativeTo"), "Q1")
    expect_equal(attr(metrics_msexp_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_msexp, "mode"), "TIC")
    expect_equal(attr(metrics_msexp_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_msexp, "change"), "jump")
    expect_equal(attr(metrics_msexp_filtered, "change"), "jump")
})
## END unit test calculateMetricsFromMsExperiment ## 

## START unit test calculateMetrics ##
## calculate the metrics by the wrapper function
suppressWarnings(
    metrics_spectra_wrapper <- calculateMetrics(object = spectra,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_spectra_wrapper_filtered <- calculateMetrics(object = spectra,
        metrics = metrics, filteredEmptySpectra = TRUE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_msexp_wrapper <- calculateMetrics(object = msexp,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))
suppressWarnings(
    metrics_msexp_wrapper_filtered <- calculateMetrics(object = msexp,
        metrics = metrics, filteredEmptySpectra = TRUE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump"))

test_that("calculateMetrics", {
    expect_equal(dim(metrics_spectra_wrapper), c(2, 12))
    expect_equal(dim(metrics_spectra_wrapper_filtered), c(2, 12))
    expect_equal(dim(metrics_msexp_wrapper), c(2, 12))
    expect_equal(dim(metrics_msexp_wrapper_filtered), c(2, 12))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_spectra_wrapper), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_spectra_wrapper_filtered), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_msexp_wrapper), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    dirs <- unlist(lapply(
        strsplit(rownames(metrics_msexp_wrapper_filtered), "sciex"), "[", 2))
    dirs <- gsub("[\\]|[/]", "", dirs)
    expect_equal(dirs, 
        c("20171016_POOL_POS_1_105-134.mzML", 
            "20171016_POOL_POS_3_105-134.mzML"))
    expect_equal(length(metrics_spectra_wrapper), 24)
    expect_equal(length(metrics_spectra_wrapper_filtered), 24)
    expect_equal(length(metrics_msexp_wrapper), 24)
    expect_equal(length(metrics_msexp_wrapper_filtered), 24)
    expect_equal(colnames(metrics_spectra_wrapper), colnames_metrics)
    expect_equal(colnames(metrics_spectra_wrapper_filtered), colnames_metrics)
    expect_equal(colnames(metrics_msexp_wrapper), colnames_metrics)
    expect_equal(colnames(metrics_msexp_wrapper_filtered), colnames_metrics)
    expect_true(is.numeric(metrics_spectra_wrapper))
    expect_true(is.numeric(metrics_spectra_wrapper_filtered))
    expect_true(is.numeric(metrics_msexp_wrapper))
    expect_true(is.numeric(metrics_msexp_wrapper_filtered))
    expect_equal(as.numeric(metrics_spectra_wrapper[1, ]), 
        metrics_spectra_1_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_wrapper_filtered[1, ]), 
        metrics_spectra_1_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_wrapper[1, ]), 
        metrics_spectra_1_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_wrapper_filtered[1, ]), 
        metrics_spectra_1_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_wrapper[2, ]), 
        metrics_spectra_2_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_wrapper_filtered[2, ]), 
        metrics_spectra_2_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_wrapper[2, ]), 
        metrics_spectra_2_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_msexp_wrapper_filtered[2, ]), 
        metrics_spectra_2_vals,  tolerance = 1e-06)
    
    expect_error(calculateMetrics(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics(spectra, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(calculateMetrics(msexp, 
        metrics = metrics, filterEmptySpectra = "foo"),
        "'filterEmptySpectra' has to be either TRUE or FALSE")
    expect_error(calculateMetrics(spectra, metrics = "foo"),
        "should be one of ")
    expect_error(calculateMetrics(msexp, metrics = "foo"),
        "should be one of ")
    
    ## test filterEmptySpectra
    expect_equal(as.numeric(calculateMetrics(object = sps_empty, 
        metrics = "numberSpectra", filterEmptySpectra = FALSE, 
        msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetrics(object = sps_empty, 
        metrics = "numberSpectra", filterEmptySpectra = TRUE, 
        msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetrics(object = sps_multiple_empty, 
        metrics = "numberSpectra", filterEmptySpectra = FALSE, 
        msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetrics(object = sps_multiple_empty, 
        metrics = "numberSpectra", filterEmptySpectra = TRUE, 
        msLevel = 2L)), 1)
    expect_equal(as.numeric(calculateMetrics(object = sps_not_empty, 
        metrics = "numberSpectra", filterEmptySpectra = FALSE, 
        msLevel = 2L)), 3)
    expect_equal(as.numeric(calculateMetrics(object = sps_not_empty, 
        metrics = "numberSpectra", filterEmptySpectra = TRUE, 
        msLevel = 2L)), 3)
    
    ## test attributes
    expect_equal(attributes(metrics_spectra_wrapper)$dimnames[[2]], 
        colnames_metrics)
    expect_equal(attributes(metrics_spectra_wrapper_filtered)$dimnames[[2]], 
        colnames_metrics)
    expect_equal(attributes(metrics_msexp_wrapper)$dimnames[[2]], 
        colnames_metrics)
    expect_equal(attributes(metrics_msexp_wrapper_filtered)$dimnames[[2]], 
        colnames_metrics)
    expect_equal(attr(metrics_spectra_wrapper, "names"), NULL)
    expect_equal(attr(metrics_spectra_wrapper_filtered, "names"), NULL)
    expect_equal(attr(metrics_msexp_wrapper, "names"), NULL)
    expect_equal(attr(metrics_msexp_wrapper_filtered, "names"), NULL)
    expect_equal(attr(metrics_spectra_wrapper, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_msexp_wrapper, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "chromatographyDuration"), 
        "MS:4000053")
    expect_equal(attr(metrics_spectra_wrapper, "ticQuartersRtFraction"),
        "MS:4000054")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_msexp_wrapper, "ticQuartersRtFraction"),
        "MS:4000054")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "ticQuartersRtFraction"), 
        "MS:4000054")
    expect_equal(attr(metrics_spectra_wrapper, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_msexp_wrapper, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "numberSpectra"), 
        "MS:4000059")
    expect_equal(attr(metrics_spectra_wrapper, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_msexp_wrapper, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "areaUnderTic"), 
        "MS:4000155")
    expect_equal(attr(metrics_spectra_wrapper, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_msexp_wrapper, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "msSignal10xChange"), 
        "MS:4000097")
    expect_equal(attr(metrics_spectra_wrapper, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_wrapper_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_msexp_wrapper, "msLevel"), 1)
    expect_equal(attr(metrics_msexp_wrapper_filtered, "msLevel"), 1)
    expect_equal(attr(metrics_spectra_wrapper, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_msexp_wrapper, "relativeTo"), "Q1")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "relativeTo"), "Q1")
    expect_equal(attr(metrics_spectra_wrapper, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_msexp_wrapper, "mode"), "TIC")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "mode"), "TIC")
    expect_equal(attr(metrics_spectra_wrapper, "change"), "jump")
    expect_equal(attr(metrics_spectra_wrapper_filtered, "change"), "jump")
    expect_equal(attr(metrics_msexp_wrapper, "change"), "jump")
    expect_equal(attr(metrics_msexp_wrapper_filtered, "change"), "jump")
})
## END unit test calculateMetrics ## 

################################################################################
################################ format = 'mzQC' ###############################
################################################################################

## calculate the metrics from Spectra
dO <- unique(spectra$dataOrigin)
spectra_1 <- spectra[spectra$dataOrigin == dO[1], ]
suppressWarnings(metrics_spectra_1 <- calculateMetricsFromOneSampleSpectra(
    spectra = spectra_1, metrics = metrics, filterEmptySpectra = FALSE, 
    msLevel = 1, relativeTo = "Q1", mode = "TIC", change = "jump", 
    format = "mzQC"))

## START unit test calculateMetricsFromOneSampleSpectra ## 
test_that("calculateMetricsFromOneSampleSpectra, format = 'mzQC'.", {
    ## spectra_1
    expect_equal(length(metrics_spectra_1), 12)
    expect_is(metrics_spectra_1, "numeric")
    ## NOTE: it should not be a list with MzQCmzQC entries
})
## END unit test calculateMetricsFromOneSampleSpectra

## START unit test calculateMetricsFromSpectra ##
## calculate the metrics from Spectra
suppressWarnings(
    metrics_spectra <- calculateMetricsFromSpectra(spectra = spectra,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump", format = "mzQC"))

test_that("calculateMetricsFromSpectra, format = 'mzQC'.", {
    expect_equal(length(metrics_spectra), 2)
    expect_equal(is(metrics_spectra[[1]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(is(metrics_spectra[[2]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(metrics_spectra[[1]]$contactAddress, as.character(NA))
    expect_equal(metrics_spectra[[2]]$contactAddress, as.character(NA))

    ## controlled vocabularies and description
    expect_equal(metrics_spectra[[1]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_spectra[[1]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_1_105-134.mzML")
    expect_equal(metrics_spectra[[2]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_spectra[[2]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_3_105-134.mzML")
    
    ## software
    expect_equal(metrics_spectra[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_spectra[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_spectra[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_spectra[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    expect_equal(metrics_spectra[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_spectra[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_spectra[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_spectra[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    
    ## chromatographyDuration
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    
    ## ticQuartersRtFraction
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    
    ## ticQuartileToQuartileLogRatio
    ## no entry for this metric since relativeTo = "Q1"
    
    ## numberSpectra
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    
    ## areaUnderTic
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.509697e+08, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.229579e+08, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    
    ## msSignal10xChange
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[1]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra[[2]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
})
## END unit test calculateMetricsFromSpectra ##


## START unit test calculateMetricsFromMsExperiment ##
msexp <- MsExperiment()
sd <- DataFrame(sample_id = c("QC1", "QC2"),
    sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
sampleData(msexp) <- sd

## define file names containing spectra data for the samples and
## add them, along with other arbitrary files to the experiment
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
experimentFiles(msexp) <- MsExperimentFiles(
    mzML_files = fls,
    annotations = "internal_standards.txt")
## link samples to data files: first sample to first file in "mzML_files",
## second sample to second file in "mzML_files"
msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
    sampleIndex = c(1, 2), withIndex = c(1, 2))
msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
    sampleIndex = c(1, 2), withIndex = c(1, 1))

## import the data and add it to the mse object
spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())

## additional parameters passed to the quality metrics functions
## (msLevel is an argument of areaUnderTic and msSignal10xChange,
## relativeTo is an argument of msSignal10xChange) passed to ...
suppressWarnings(metrics_msexp <- calculateMetricsFromMsExperiment(msexp = msexp, 
    metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
    relativeTo = "Q1", mode = "TIC", change = "jump", format = "mzQC"))

test_that("calculateMetricsFromMsExperiment, format = 'mzQC'.", {
    expect_equal(length(metrics_msexp), 2)
    expect_equal(is(metrics_msexp[[1]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(is(metrics_msexp[[2]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(metrics_msexp[[1]]$contactAddress, as.character(NA))
    expect_equal(metrics_msexp[[2]]$contactAddress, as.character(NA))
    
    
    ## controlled vocabularies and description
    expect_equal(metrics_msexp[[1]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_msexp[[1]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_1_105-134.mzML")
    expect_equal(metrics_msexp[[2]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_msexp[[2]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_3_105-134.mzML")
    
    ## software
    expect_equal(metrics_msexp[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_msexp[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_msexp[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_msexp[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    expect_equal(metrics_msexp[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_msexp[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_msexp[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_msexp[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    
    ## chromatographyDuration
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    
    ## ticQuartersRtFraction
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    
    ## ticQuartileToQuartileLogRatio
    ## no entry for this metric since relativeTo = "Q1"
    
    ## numberSpectra
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    
    ## areaUnderTic
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.509697e+08, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.229579e+08, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    
    ## msSignal10xChange
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[1]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp[[2]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
})
## END unit test calculateMetricsFromMsExperiment ## 

## START unit test calculateMetrics ##
## calculate the metrics by the wrapper function
suppressWarnings(
    metrics_spectra_wrapper <- calculateMetrics(object = spectra,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump", format = "mzQC"))
suppressWarnings(
    metrics_msexp_wrapper <- calculateMetrics(object = msexp,
        metrics = metrics, filterEmptySpectra = FALSE, msLevel = 1, 
        relativeTo = "Q1", mode = "TIC", change = "jump", format = "mzQC"))

test_that("calculateMetrics, format = 'mzQC'.", {
    
    ## metrics_spectra_wrapper
    expect_equal(length(metrics_spectra_wrapper), 2)
    expect_equal(is(metrics_spectra_wrapper[[1]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(is(metrics_spectra_wrapper[[2]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(metrics_spectra_wrapper[[1]]$contactAddress, as.character(NA))
    expect_equal(metrics_spectra_wrapper[[2]]$contactAddress, as.character(NA))
    
    ## controlled vocabularies and description
    expect_equal(metrics_spectra_wrapper[[1]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_spectra_wrapper[[1]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_1_105-134.mzML")
    expect_equal(metrics_spectra_wrapper[[2]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_spectra_wrapper[[2]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_3_105-134.mzML")
    
    ## software
    expect_equal(metrics_spectra_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_spectra_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_spectra_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_spectra_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
         "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    expect_equal(metrics_spectra_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_spectra_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_spectra_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_spectra_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    
    ## chromatographyDuration
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    
    ## ticQuartersRtFraction
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    
    ## ticQuartileToQuartileLogRatio
    ## no entry for this metric since relativeTo = "Q1"
    
    ## numberSpectra
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    
    ## areaUnderTic
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.509697e+08, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.229579e+08, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    
    ## msSignal10xChange
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_spectra_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
    
    ## 
    ## metrics_msexp_wrapper
    expect_equal(length(metrics_msexp_wrapper), 2)
    expect_equal(is(metrics_msexp_wrapper[[1]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(is(metrics_msexp_wrapper[[2]]), c("MzQCmzQC", "envRefClass", 
        ".environment", "refClass", "environment", "refObject"))
    expect_equal(metrics_msexp_wrapper[[1]]$contactAddress, as.character(NA))
    expect_equal(metrics_msexp_wrapper[[2]]$contactAddress, as.character(NA))
    
    ## controlled vocabularies and description
    expect_equal(metrics_msexp_wrapper[[1]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_msexp_wrapper[[1]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_1_105-134.mzML")
    expect_equal(metrics_msexp_wrapper[[2]]$controlledVocabularies[[1]]$name, 
        "Proteomics Standards Initiative Mass Spectrometry Ontology")
    expect_equal(metrics_msexp_wrapper[[2]]$description,
        "A mzQC document on the sample 20171016_POOL_POS_3_105-134.mzML")
    
    ## software
    expect_equal(metrics_msexp_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_msexp_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_msexp_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_msexp_wrapper[[1]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    expect_equal(metrics_msexp_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$accession, "MS:4000151")
    expect_equal(metrics_msexp_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$name, "MsQuality")
    expect_equal(metrics_msexp_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$version, 
        packageDescription("MsQuality")$Version)
    expect_equal(metrics_msexp_wrapper[[2]]$runQualities[[1]]$metadata$analysisSoftware[[1]]$description, 
        "\"MsQuality – an interoperable open-source package for the calculation of standardized quality metrics of mass spectrometry data.\" [DOI:10.1101/2023.05.12.540477, https://github.com/tnaake/MsQuality/]")
    
    ## chromatographyDuration
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$accession, 
        "MS:4000053")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$name, 
        "chromatography duration")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$description, 
        "\"The retention time duration of the chromatography in seconds.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$value, 
        2.594770e+02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[1]]$unit, 
        list())
    
    ## ticQuartersRtFraction
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$accession, 
        "MS:4000054")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$name, 
        "TIC quarters RT fraction")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$description, 
        "\"The interval when the respective quarter of the TIC accumulates divided by retention time duration.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[2]]$unit, 
        list())
    
    ## ticQuartileToQuartileLogRatio
    ## no entry for this metric since relativeTo = "Q1"
    
    ## numberSpectra
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$accession, 
        "MS:4000059")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$name, 
        "number of MS1 spectra")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$description, 
        "\"The number of MS1 events in the run.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$value, 
        9.310000e02, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[3]]$unit, 
        list())
    
    ## areaUnderTic
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.509697e+08, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$accession,
        "MS:4000155")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$name,
        "area under TIC")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$description,
        "\"The area under the total ion chromatogram.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$value,
        6.229579e+08, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[4]]$unit,
        list())
    
    ## msSignal10xChange
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[1]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$accession, 
        "MS:4000097")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$name, 
        "MS1 signal jump (10x) count")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$description, 
        "\"The number of times where MS1 TIC increased more than 10-fold between adjacent MS1 scans.\" [PSI:MS]")
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$value, 
        0, tolerance = 1e-06)
    expect_equal(
        metrics_msexp_wrapper[[2]]$runQualities[[1]]$qualityMetrics[[5]]$unit, 
        list())
})
## END unit test calculateMetrics ## 

