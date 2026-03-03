## create toy example data set (Spectra)
library(MsDataHub)
library(Spectra)
library(Chromatograms)
library(MsExperiment)
library(S4Vectors)
## define file names containing spectra data for the samples
sciex_file <- c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
    MsDataHub::X20171016_POOL_POS_3_105.134.mzML())

## import the data and assign it to the spectra object
spectra <- Spectra(files = sciex_file)
chr <- Chromatograms(spectra, "sum")
## create toy example data set (MsExperiment)
msexp <- MsExperiment()
sd <- DataFrame(sample_id = c("QC1", "QC2"),
                sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
sampleData(msexp) <- sd

## define file names containing spectra data for the samples and
## add them, along with other arbitrary files to the experiment
experimentFiles(msexp) <- MsExperimentFiles(
    mzML_files = sciex_file,
    annotations = "internal_standards.txt")
## link samples to data files: first sample to first file in "mzML_files",
## second sample to second file in "mzML_files"
msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
    sampleIndex = c(1, 2), withIndex = c(1, 2))
msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
    sampleIndex = c(1, 2), withIndex = c(1, 1))

## import the data and add it to the msexp object
spectra(msexp) <- Spectra(sciex_file)


qm_spectra <- c("chromatographyDuration", "ticQuantileRtFraction",
    "rtOverMsQuarters", "ticQuartileToQuartileLogRatio", "numberSpectra",
    "numberEmptyScans", "medianPrecursorMz", "rtIqr", "rtIqrRate",
    "areaUnderTic", "areaUnderTicRtQuantiles",
    "extentIdentifiedPrecursorIntensity", "medianTicRtIqr",
    "medianTicOfRtRange", "mzAcquisitionRange", "rtAcquisitionRange",
    "precursorIntensityRange", "precursorIntensityQuartiles",
    "precursorIntensityMean", "precursorIntensitySd", "msSignal10xChange",
    "ratioCharge1over2", "ratioCharge3over2", "ratioCharge4over2", "meanCharge",
    "medianCharge")
qm_mse <- c("chromatographyDuration", "ticQuantileRtFraction",
    "rtOverMsQuarters", "ticQuartileToQuartileLogRatio", "numberSpectra",
    "numberEmptyScans", "medianPrecursorMz", "rtIqr", "rtIqrRate",
    "areaUnderTic", "areaUnderTicRtQuantiles",
    "extentIdentifiedPrecursorIntensity", "medianTicRtIqr",
    "medianTicOfRtRange", "mzAcquisitionRange", "rtAcquisitionRange",
    "precursorIntensityRange", "precursorIntensityQuartiles",
    "precursorIntensityMean", "precursorIntensitySd", "msSignal10xChange",
    "ratioCharge1over2", "ratioCharge3over2", "ratioCharge4over2",
    "meanCharge", "medianCharge")
qm_chr <-  c("chromatographyDuration", "chromatogramCount",
            "rtAcquisitionRange", "maxIntensity",
            "intensityQuartiles", "intensityMean", "intensitySd",
            "intensityRange", "peakCount", "rtIqr",
            "baselineIntensity", "signalToNoiseRatio",
            "intensity10xChange",
            "numberEmptyChrom",
            "medianIntensityRtIqr",
            "xicFwhm", "peakBoundary", "peakWidth", "peakBeta",
            "peakProminence",
            "ticQuantileRtFraction", "areaUnderTic", "areaUnderTicRtQuantiles",
            "areaUnderTicMs1", "areaUnderTicMs2")

test_that("qualityMetrics for Spectra", {
    expect_equal(qualityMetrics(spectra), qm_spectra)
    expect_equal(qualityMetrics(spectra), qm_mse)
    expect_error(qualityMetrics(NULL), "object '.metrics' not found")
})

test_that("qualityMetrics for Chromatograms", {
    chr <- Chromatograms(spectra)
    expect_equal(qualityMetrics(chr), qm_chr)
})

test_that("qualityMetrics for MsExperiment", {
    expect_equal(qualityMetrics(msexp), qm_mse)
    expect_is(msexp, "MsExperiment")
    expect_equal(length(spectra(msexp)), 1862)
})
