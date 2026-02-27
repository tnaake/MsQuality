## create toy example data set (Spectra)
library(MsDataHub)
library(Spectra)
## define file names containing spectra data for the samples
sciex_file <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
    MsDataHub::X20171016_POOL_POS_3_105.134.mzML())

## import the data and assign it to the spectra object
spectra <- Spectra(sciex_file)

## define the quality metrics to be calculated
metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")

## calculate the metrics
## additional parameters passed to the quality metrics functions
## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
## relativeTo is an argument of msSignal10XChange)
qc <- calculateMetricsFromSpectra(spectra = spectra, metrics = metrics,
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
    expect_equal(tbl$value, c(650969677, 622957884))
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
