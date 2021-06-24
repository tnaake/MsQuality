#' @name qualityMetrics
#'
#' @title Get a vector of quality metrics than can be applied to `object`
#'
#' @description
#' The function `qualityMetrics` returns a character vector with available
#' quality metrics depending on `object`.
#' 
#' @details
#' `object` is either a `Spectra` or `MsExperiment` object. 
#' 
#' @param object object of type `Spectra` or `MsExperiment`
#' 
#' @return `character`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @examples 
#' library(msdata)
#' library(Spectra)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' qualityMetrics(spectra)
qualityMetrics <- function(object) {
    if (is(object, "Spectra"))
        .metrics <- c(
            "rtDuration", "rtOverTICquantile", "rtOverMSQuarters",
            "ticQuantileToQuantileLogRatio", "numberSpectra", 
            "medianPrecursorMZ", "rtIQR", "rtIQRrate", "areaUnderTIC", 
            "areaUnderTICRTquantiles", "extentIdentifiedPrecursorIntensity",
            "medianTICRTIQR", "medianTICofRTRange", "mzAcquisitionRange",
            "rtAcquisitionRange", "precursorIntensityRange", 
            "precursorIntensityQuartiles", "precursorIntensityMean",
            "precursorIntensitySD", "msSignal10XChange", "RatioCharge1over2", 
            "RatioCharge3over2", "RatioCharge4over2", "meanCharge", 
            "medianCharge"
        )
    if (is(object, "MsExperiment"))
        .metrics <-  c(
            "rtDuration", "rtOverTICquantile", "rtOverMSQuarters",
            "ticQuantileToQuantileLogRatio", "numberSpectra", 
            "medianPrecursorMZ", "rtIQR", "rtIQRrate", "areaUnderTIC", 
            "areaUnderTICRTquantiles", "extentIdentifiedPrecursorIntensity", 
            "medianTICRTIQR", "medianTICofRTRange", "mzAcquisitionRange",
            "rtAcquisitionRange", "precursorIntensityRange", 
            "precursorIntensityQuartiles", "precursorIntensityMean",
            "precursorIntensitySD", "msSignal10XChange", "RatioCharge1over2", 
            "RatioCharge3over2", "RatioCharge4over2", "meanCharge", 
            "medianCharge"
        )
    return(.metrics)
}

