## create toy example data set (Spectra)
library(msdata)
library(Spectra)
fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
spectra <- Spectra(fls, backend = MsBackendMzR())

qm_spectra <- c("rtDuration", "rtOverTicQuantile", "rtOverMsQuarters",
    "ticQuartileToQuartileLogRatio", "numberSpectra", "medianPrecursorMz", 
    "rtIqr", "rtIqrRate", "areaUnderTic", "areaUnderTicRtQuantiles", 
    "extentIdentifiedPrecursorIntensity", "medianTicRtIqr", 
    "medianTicOfRtRange", "mzAcquisitionRange", "rtAcquisitionRange", 
    "precursorIntensityRange", "precursorIntensityQuartiles", 
    "precursorIntensityMean", "precursorIntensitySd", "msSignal10xChange", 
    "ratioCharge1over2", "ratioCharge3over2", "ratioCharge4over2", "meanCharge", 
    "medianCharge")

test_that("qualityMetrics", {
    expect_equal(qualityMetrics(spectra), qm_spectra)
    expect_error(qualityMetrics(NULL), "object '.metrics' not found")
})

