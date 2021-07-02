## create Spectra object
library(msdata)
library(Spectra)
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## create MsExperiment object
library(msdata)
library(MsExperiment)
library(S4Vectors)
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
mse <- linkSampleData(mse, with = "experimentFiles.annotations",
    sampleIndex = c(1, 2), withIndex = c(1, 1))

## import the data and add it to the mse object
spectra(mse) <- Spectra(fls, backend = MsBackendMzR())
 
## build the results
## define the quality metrics to be calculated
metrics <- c("rtDuration", "rtOverTICquantile", "ticQuantileToQuantileLogRatio", 
    "numberSpectra", "areaUnderTIC", "msSignal10Change")

## additional parameters passed to the quality metrics functions
## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
    change = c("jump", "fall"))

## calculate the metrics from MsExperiment
metrics_mse <- calculateMetricsFromMsExperiment(mse = mse, 
    metrics = metrics, params = params_l)

## calculate the metrics from Spectra
metrics_spectra <- calculateMetricsFromSpectra(spectra = spectra, 
    metrics = metrics, params = params_l)

## calculate the metrics by the wrapper function
metrics_spectra_wrapper <- calculateMetrics(object = spectra, 
    metrics = metrics, params = params_l)
metrics_mse_wrapper <- calculateMetrics(object = mse, 
    metrics = metrics, params = params_l)

## START unit test calculateMetricsFromMsExperiment ## 
colnames_metrics_mse <- c("rtDuration", "rtOverTICquantile_MSLevel1_0%",                                  
    "rtOverTICquantile_MSLevel1_25%", "rtOverTICquantile_MSLevel1_50%",                               
    "rtOverTICquantile_MSLevel1_75%", "rtOverTICquantile_MSLevel1_100%",                             
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToQ1_Q2/Q1", 
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToQ1_Q3/Q1",
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToQ1_Q4/Q1", 
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToprevious_Q2/Q1",
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToprevious_Q3/Q2", 
    "ticQuantileToQuantileLogRatio_MSLevel1_relativeToprevious_Q4/Q3",
    "numberSpectra_MSLevel1", "areaUnderTIC_MSLevel1")
test_that("calculateMetricsFromMsExperiment", {
    expect_equal(dim(metrics_mse), c(2, 14))
    expect_equal(rownames(metrics_mse), NULL)
    expect_equal(colnames(metrics_mse), colnames_metrics_mse)
    expect_true(is.numeric(metrics_mse))
    expect_error(calculateMetricsFromMsExperiment(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromMsExperiment("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromMsExperiment(mse, metrics = "foo"),
        "should be one of")
})
## END unit test calculateMetricsFromMsExperiment

## START unit test calculateMetricsFromSpectra ##
metrics_spectra_vals <- c(2.594820e+02, 1.058682e-03, 2.502724e-01, 
    5.005178e-01, 7.497315e-01, 1.000000e+00, -5.853477e-02, -6.405647e-01, 
    -4.869522e-01, -5.853477e-02, -5.820299e-01, 1.536125e-01, 1.862000e+03,
    1.273928e+09)
test_that("calculateMetricsFromSpectra", {
    expect_equal(length(metrics_spectra), 14)
    expect_equal(names(metrics_spectra), colnames_metrics_mse)
    expect_true(is.numeric(metrics_spectra))
    expect_equal(as.numeric(metrics_spectra), metrics_spectra_vals, 
        tolerance = 1e-06)
    expect_error(calculateMetricsFromSpectra(NULL, metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromSpectra("foo", metrics = metrics),
        "object '.metrics' not found")
    expect_error(calculateMetricsFromSpectra(spectra, metrics = "foo"),
        "should be one of ")
})
## END unit test calculateMetricsFromSpectra ##

## START unit test calculateMetrics ##
test_that("calculateMetrics", {
    expect_equal(length(metrics_spectra_wrapper), 14)
    expect_equal(names(metrics_spectra_wrapper), colnames_metrics_mse)
    expect_true(is.numeric(metrics_spectra_wrapper))
    expect_equal(as.numeric(metrics_spectra_wrapper), metrics_spectra_vals,
        tolerance = 1e-06)
    expect_equal(dim(metrics_mse_wrapper), c(2, 14))
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
