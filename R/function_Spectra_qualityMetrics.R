#' @name qualityMetrics
#'
#' @title Get a vector of quality metrics than can be applied to \code{object}
#'
#' @description
#' The function \code{qualityMetrics} returns a character vector with available
#' quality metrics depending on \code{object}.
#' 
#' @details
#' \code{object} is a \code{Spectra} or \code{MsExperiment}. 
#' 
#' @param object object of type \code{Spectra} or \code{MsExperiment}
#' 
#' @return \code{character}
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom Spectra Spectra
#' @importFrom MsExperiment MsExperiment
#' 
#' @export
#' 
#' @examples 
#' library(Spectra)
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$dataOrigin <- rep("sample_1", 3)
#' sps <- Spectra(spd)
#' 
#' qualityMetrics(object = sps)
qualityMetrics <- function(object) {
    if (is(object, "Spectra"))
        .metrics <- c(
            "rtDuration", "rtOverTicQuantiles", "rtOverMsQuarters",
            "ticQuartileToQuartileLogRatio", "numberSpectra", 
            "medianPrecursorMz", "rtIqr", "rtIqrRate", "areaUnderTic", 
            "areaUnderTicRtQuantiles", "extentIdentifiedPrecursorIntensity",
            "medianTicRtIqr", "medianTicOfRtRange", "mzAcquisitionRange",
            "rtAcquisitionRange", "precursorIntensityRange", 
            "precursorIntensityQuartiles", "precursorIntensityMean",
            "precursorIntensitySd", "msSignal10xChange", "ratioCharge1over2", 
            "ratioCharge3over2", "ratioCharge4over2", "meanCharge", 
            "medianCharge"
        )
    
    if (is(object, "MsExperiment"))
        .metrics <-  c(
            "rtDuration", "rtOverTicQuantiles", "rtOverMsQuarters",
            "ticQuartileToQuartileLogRatio", "numberSpectra", 
            "medianPrecursorMz", "rtIqr", "rtIqrRate", "areaUnderTic", 
            "areaUnderTicRtQuantiles", "extentIdentifiedPrecursorIntensity", 
            "medianTicRtIqr", "medianTicOfRtRange", "mzAcquisitionRange",
            "rtAcquisitionRange", "precursorIntensityRange", 
            "precursorIntensityQuartiles", "precursorIntensityMean",
            "precursorIntensitySd", "msSignal10xChange", "ratioCharge1over2", 
            "ratioCharge3over2", "ratioCharge4over2", "meanCharge", 
            "medianCharge"
      )
    
    .metrics
}
