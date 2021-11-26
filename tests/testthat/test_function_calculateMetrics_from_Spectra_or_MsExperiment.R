fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## create MsExperiment object
mse <- MsExperiment()
sd <- DataFrame(sample_id = c("QC1", "QC2"),
                sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
sampleData(mse) <- sd

## define file names containing spectra data for the samples and
## add them, along with other arbitrary files to the experiment
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
experimentFiles(mse) <- MsExperimentFiles(
    mzML_files = fls,
    annotations = "internal_standards.txt")
## link samples to data files: first sample to first file in "mzML_files",
## second sample to second file in "mzML_files"
mse <- linkSampleData(mse, with = "experimentFiles.mzML_files",
                      sampleIndex = c(1, 2), withIndex = c(1, 2))

## import the data and add it to the mse object
spectra(mse) <- spectra

## build the results
## define the quality metrics to be calculated
metrics <- c("rtDuration", "rtOverTicQuantile", "ticQuantileToQuantileLogRatio",
    "numberSpectra", "areaUnderTic", "msSignal10xChange")

## calculate the metrics from MsExperiment
## additional parameters passed to the quality metrics functions
## (msLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
suppressWarnings(metrics_mse <- calculateMetricsFromMsExperiment(mse = mse,
    metrics = metrics, msLevel = 1, relativeTo = "Q1", change = "jump"))

## calculate the metrics from Spectra
suppressWarnings(metrics_spectra <- calculateMetricsFromSpectra(
    spectra = spectra, metrics = metrics, msLevel = 1, 
    relativeTo = "Q1", change = "jump"))

## calculate the metrics by the wrapper function
suppressWarnings(metrics_spectra_wrapper <- calculateMetrics(object = spectra,
    metrics = metrics, msLevel = 1, relativeTo = "Q1", change = "jump"))
suppressWarnings(metrics_mse_wrapper <- calculateMetrics(object = mse,
    metrics = metrics, msLevel = 1, relativeTo = "Q1", change = "jump"))

## START unit test calculateMetricsFromMsExperiment ## 
colnames_metrics_mse <- c("rtDuration", "rtOverTicQuantile.0%",                  
    "rtOverTicQuantile.25%", "rtOverTicQuantile.50%",
    "rtOverTicQuantile.75%", "rtOverTicQuantile.100%",
    "ticQuantileToQuantileLogRatio.Q2/Q1",
    "ticQuantileToQuantileLogRatio.Q3/Q1",
    "ticQuantileToQuantileLogRatio.Q4/Q1",
    "numberSpectra", "areaUnderTic", "msSignal10xChange")
test_that("calculateMetricsFromMsExperiment", {
    expect_equal(dim(metrics_mse), c(2, 12))
    expect_equal(rownames(metrics_mse), NULL)
    expect_equal(colnames(metrics_mse), colnames_metrics_mse)
    expect_true(is.numeric(metrics_mse))
    expect_error(calculateMetricsFromMsExperiment(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromMsExperiment("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromMsExperiment(mse, metrics = "foo"),
        "should be one of")
    expect_error(calculateMetricsFromMsExperiment(mse, 
        metrics = "ticQuantileToQuantileLogRatio", 
        relativeTo = c("Q1", "previous")), "'relativeTo' has to be of length 1")
    expect_error(calculateMetricsFromMsExperiment(mse, 
        metrics = "msSignal10xChange", change = c("jump", "fall")), 
        "'change' has to be of length 1")
})
## END unit test calculateMetricsFromMsExperiment

## START unit test calculateMetricsFromSpectra ##
metrics_spectra_vals <- c(2.594820e+02, 0, 2.505338e-01,
    5.000077e-01, 7.505222e-01, 1.000000e+00, -5.853477e-02, -6.405647e-01,
    -4.869522e-01, 1.862000e+03, 1.273928e+09, 0)
test_that("calculateMetricsFromSpectra", {
    expect_equal(length(metrics_spectra), 12)
    expect_equal(names(metrics_spectra), colnames_metrics_mse)
    expect_true(is.numeric(metrics_spectra))
    expect_equal(as.numeric(metrics_spectra), metrics_spectra_vals, 
        tolerance = 1e-06)
    expect_equal(attributes(metrics_spectra)$names, colnames_metrics_mse)
    expect_equal(attributes(metrics_spectra)$msLevel, 1)
    expect_equal(attributes(metrics_spectra)$relativeTo, "Q1")
    expect_equal(attributes(metrics_spectra)$change, "jump")
    expect_error(calculateMetricsFromSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
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
    expect_equal(length(metrics_spectra_wrapper), 12)
    expect_equal(names(metrics_spectra_wrapper), colnames_metrics_mse)
    expect_true(is.numeric(metrics_spectra_wrapper))
    expect_equal(as.numeric(metrics_spectra_wrapper), metrics_spectra_vals,
        tolerance = 1e-06)
    expect_equal(dim(metrics_mse_wrapper), c(2, 12))
    expect_equal(rownames(metrics_mse_wrapper), NULL)
    expect_equal(colnames(metrics_mse_wrapper), colnames_metrics_mse)
    expect_true(is.numeric(metrics_mse_wrapper))
    expect_error(calculateMetrics(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetrics(spectra, metrics = "foo"),
        "should be one of ")
    expect_error(calculateMetrics(mse, metrics = "foo"),
        "should be one of ")
    })
## END unit test calculateMetrics ## 

