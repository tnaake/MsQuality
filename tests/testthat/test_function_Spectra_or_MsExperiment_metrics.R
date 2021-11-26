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


qm_spectra <- c("rtDuration", "rtOverTicQuantile", "rtOverMsQuarters",
    "ticQuantileToQuantileLogRatio", "numberSpectra", "medianPrecursorMz", 
    "rtIqr", "rtIqrRate", "areaUnderTic", "areaUnderTicRtQuantiles", 
    "extentIdentifiedPrecursorIntensity", "medianTicRtIqr", 
    "medianTicOfRtRange", "mzAcquisitionRange", "rtAcquisitionRange", 
    "precursorIntensityRange", "precursorIntensityQuartiles", 
    "precursorIntensityMean", "precursorIntensitySd", "msSignal10xChange", 
    "ratioCharge1over2", "ratioCharge3over2", "ratioCharge4over2", "meanCharge", 
    "medianCharge")
qm_mse <- c("rtDuration", "rtOverTicQuantile", "rtOverMsQuarters",
    "ticQuantileToQuantileLogRatio", "numberSpectra", "medianPrecursorMz", 
    "rtIqr", "rtIqrRate", "areaUnderTic", "areaUnderTicRtQuantiles", 
    "extentIdentifiedPrecursorIntensity", "medianTicRtIqr", 
    "medianTicOfRtRange", "mzAcquisitionRange", "rtAcquisitionRange", 
    "precursorIntensityRange", "precursorIntensityQuartiles", 
    "precursorIntensityMean", "precursorIntensitySd", "msSignal10xChange", 
    "ratioCharge1over2", "ratioCharge3over2", "ratioCharge4over2", "meanCharge", 
    "medianCharge")

test_that("qualityMetrics", {
    expect_equal(qualityMetrics(spectra), qm_spectra)
    expect_equal(qualityMetrics(spectra), qm_mse)
    expect_error(qualityMetrics(NULL), "object '.metrics' not found")
})

