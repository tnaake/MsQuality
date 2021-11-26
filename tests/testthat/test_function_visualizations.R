## create toy example data set (Spectra)
library(msdata)
library(Spectra)
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

## create toy example data set (MsExperiment)
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

library(Spectra)
## import the data and add it to the mse object
spectra(mse) <- Spectra(fls, backend = MsBackendMzR())

## define the quality metrics to be calculated
metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")

## calculate the metrics
## additional parameters passed to the quality metrics functions
## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
qc <- calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
    msLevel = 1, relativeTo = "Q1", change = "jump")
rownames(qc) <- c("Sample 1", "Sample 2")

## START unit test plotMetric ## 
test_that("plotMetric", {
    expect_is(plotMetric(qc = qc, metric = "areaUnderTic"), "plotly")
    expect_error(plotMetric(qc = NULL, metric = "areaUnderTic"),
        "'metric' not in qc")
    expect_error(plotMetric(qc = matrix(), metric = "areaUnderTic"),
        "'metric' not in qc")
    expect_error(plotMetric(qc = qc, metric = "foo"), "'metric' not in qc")
})
## END unit test plotMetric ##

## START unit test plotMetricTibble ##
test_that("plotMetricTibble", {
    expect_is(plotMetricTibble(qc = qc, metric = "areaUnderTic"), "tbl")
    tbl <- plotMetricTibble(qc = qc, metric = "areaUnderTic")
    expect_equal(tbl$rowname, factor(c("Sample 1", "Sample 2")))
    expect_equal(tbl$name, c("areaUnderTic", "areaUnderTic"))
    expect_equal(tbl$value, c(1273927561, 1273927561))
    expect_error(plotMetricTibble(qc = NULL, metric = "areaUnderTic"),
                 "'metric' not in qc")
    expect_error(plotMetricTibble(qc = matrix(), metric = "areaUnderTic"),
                 "'metric' not in qc")
    expect_error(plotMetricTibble(qc = qc, metric = "foo"), "'metric' not in qc")
})
## END unit test plotMetricTibble ## 

## START unit test shinyMsQuality ##
test_that("shinyMsQuality", {
    expect_error(shinyMsQuality(qc = matrix()), "'qc' has to be numeric")
    expect_error(shinyMsQuality(qc = NULL), "'qc' is not a matrix")
})
## END unit test shinyMsQuality ## 
