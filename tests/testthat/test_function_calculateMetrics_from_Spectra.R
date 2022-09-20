fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## build the results
## define the quality metrics to be calculated
metrics <- c("rtDuration", "rtOverTicQuantile", "ticQuantileToQuantileLogRatio",
    "numberSpectra", "areaUnderTic", "msSignal10xChange")

## additional parameters passed to the quality metrics functions
## (msLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
suppressWarnings(metrics_spectra <- calculateMetricsFromSpectra(spectra = spectra,
    metrics = metrics, msLevel = 1, relativeTo = "Q1", change = "jump"))

## calculate the metrics from Spectra
dO <- unique(spectra$dataOrigin)
spectra_1 <- spectra[spectra$dataOrigin == dO[1], ]
spectra_2 <- spectra[spectra$dataOrigin == dO[2], ]
suppressWarnings(metrics_spectra_1 <- calculateMetricsFromOneSampleSpectra(
    spectra = spectra_1, metrics = metrics, msLevel = 1, 
    relativeTo = "Q1", change = "jump"))
suppressWarnings(metrics_spectra_2 <- calculateMetricsFromOneSampleSpectra(
    spectra = spectra_2, metrics = metrics, msLevel = 1, 
    relativeTo = "Q1", change = "jump"))

## calculate the metrics by the wrapper function
suppressWarnings(metrics_spectra_wrapper <- calculateMetrics(object = spectra,
    metrics = metrics, msLevel = 1, relativeTo = "Q1", change = "jump"))

## START unit test calculateMetricsFromOneSampleSpectra ## 
colnames_metrics <- c("rtDuration", "rtOverTicQuantile.0%",                  
    "rtOverTicQuantile.25%", "rtOverTicQuantile.50%",
    "rtOverTicQuantile.75%", "rtOverTicQuantile.100%",
    "ticQuantileToQuantileLogRatio.Q2/Q1",
    "ticQuantileToQuantileLogRatio.Q3/Q1",
    "ticQuantileToQuantileLogRatio.Q4/Q1",
    "numberSpectra", "areaUnderTic", "msSignal10xChange")
metrics_spectra_1_vals <- c(2.594770e+02, 0, 2.505386e-01,
    4.999981e-01, 7.505367e-01, 1.000000e+00, -1.603689e-01, -7.040791e-01,
    -5.018226e-01, 9.310000e02, 6.509697e+08, 0)
metrics_spectra_2_vals <- c(2.594770e+02, 0, 2.505386e-01,
    4.999981e-01, 7.505367e-01, 1.000000e+00, 5.052683e-02, -5.673914e-01,
    -4.700149e-01, 9.310000e02, 6.229579e+08, 0)

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
    expect_equal(attributes(metrics_spectra_1)$names, colnames_metrics)
    expect_equal(attributes(metrics_spectra_1)$msLevel, 1)
    expect_equal(attributes(metrics_spectra_1)$relativeTo, "Q1")
    expect_equal(attributes(metrics_spectra_1)$change, "jump")
    expect_error(calculateMetricsFromOneSampleSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromOneSampleSpectra("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(
        calculateMetricsFromOneSampleSpectra(spectra_1, metrics = "foo"),
        "should be one of")
    expect_error(
        calculateMetricsFromOneSampleSpectra(spectra_1, 
        metrics = "ticQuantileToQuantileLogRatio", 
        relativeTo = c("Q1", "previous")), "'relativeTo' has to be of length 1")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_1, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
    
    ## spectra_2
    expect_true(is.numeric(metrics_spectra_2))
    expect_equal(length(metrics_spectra_2), 12)
    expect_equal(names(metrics_spectra_2), colnames_metrics)
    expect_true(is.numeric(metrics_spectra_2))
    expect_equal(as.numeric(metrics_spectra_2), metrics_spectra_2_vals, 
        tolerance = 1e-06)
    expect_equal(attributes(metrics_spectra_2)$names, colnames_metrics)
    expect_equal(attributes(metrics_spectra_2)$msLevel, 1)
    expect_equal(attributes(metrics_spectra_2)$relativeTo, "Q1")
    expect_equal(attributes(metrics_spectra_2)$change, "jump")
    expect_error(calculateMetricsFromOneSampleSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_2,
        metrics = "ticQuantileToQuantileLogRatio",
        relativeTo = c("Q1", "previous")), "'relativeTo' has to be of length 1")
    expect_error(calculateMetricsFromOneSampleSpectra(spectra_2,
        metrics = "msSignal10xChange", change = c("jump", "fall")),
        "'change' has to be of length 1")
})
## END unit test calculateMetricsFromOneSampleSpectra

## START unit test calculateMetricsFromSpectra ##
test_that("calculateMetricsFromSpectra", {
    expect_equal(dim(metrics_spectra), c(2, 12))
    expect_equal(unlist(lapply(
            strsplit(rownames(metrics_spectra), "sciex"), "[", 2)), 
        c("\\20171016_POOL_POS_1_105-134.mzML", 
          "\\20171016_POOL_POS_3_105-134.mzML"))
    expect_equal(colnames(metrics_spectra), colnames_metrics)
    expect_equal(as.numeric(metrics_spectra[1, ]), 
        metrics_spectra_1_vals, tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra[2, ]), 
        metrics_spectra_2_vals, tolerance = 1e-06)
    expect_equal(attributes(metrics_spectra)$dimnames[[2]], colnames_metrics)
    expect_equal(attributes(metrics_spectra)$msLevel, 1)
    expect_equal(attributes(metrics_spectra)$relativeTo, "Q1")
    expect_equal(attributes(metrics_spectra)$change, "jump")
    
    expect_error(calculateMetricsFromSpectra("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromSpectra(spectra, metrics = "foo"),
        "should be one of ")
    expect_error(calculateMetricsFromSpectra(spectra, 
        metrics = "ticQuantileToQuantileLogRatio",
        relativeTo = c("Q1", "previous")), "'relativeTo' has to be of length 1")
    expect_error(calculateMetricsFromSpectra(spectra, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
    
})
## END unit test calculateMetricsFromSpectra ##

## START unit test calculateMetrics ##
test_that("calculateMetrics", {
    expect_equal(dim(metrics_spectra_wrapper), c(2, 12))
    expect_equal(unlist(lapply(
            strsplit(rownames(metrics_spectra_wrapper), "sciex"), "[", 2)), 
        c("\\20171016_POOL_POS_1_105-134.mzML", 
            "\\20171016_POOL_POS_3_105-134.mzML"))
    expect_equal(length(metrics_spectra_wrapper), 24)
    expect_equal(colnames(metrics_spectra_wrapper), colnames_metrics)
    expect_true(is.numeric(metrics_spectra_wrapper))
    expect_equal(as.numeric(metrics_spectra_wrapper[1, ]), 
        metrics_spectra_1_vals,  tolerance = 1e-06)
    expect_equal(as.numeric(metrics_spectra_wrapper[2, ]), 
        metrics_spectra_2_vals,  tolerance = 1e-06)
    expect_equal(attributes(metrics_spectra_wrapper)$dimnames[[2]], 
        colnames_metrics)
    expect_equal(attributes(metrics_spectra_wrapper)$msLevel, 1)
    expect_equal(attributes(metrics_spectra_wrapper)$relativeTo, "Q1")
    expect_equal(attributes(metrics_spectra_wrapper)$change, "jump")
    
    expect_error(calculateMetrics(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics(spectra, metrics = "foo"),
        "should be one of ")
})
## END unit test calculateMetrics ## 

