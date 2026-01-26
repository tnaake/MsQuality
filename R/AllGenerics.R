## S4 Generic definitions for MsQuality
## These must be defined before their methods in other files

#' @rdname calculateMetrics
#' @export
setGeneric("calculateMetrics", function(object, metrics, ...) {
    standardGeneric("calculateMetrics")
})

#' @rdname areaUnderTic
#' @export
setGeneric("areaUnderTic", function(object, ...) standardGeneric("areaUnderTic"))

#' @rdname ticQuantileRtFraction
#' @export
setGeneric("ticQuantileRtFraction", function(object, ...) standardGeneric("ticQuantileRtFraction"))

#' @rdname areaUnderTicRtQuantiles
#' @export
setGeneric("areaUnderTicRtQuantiles", function(object, ...) standardGeneric("areaUnderTicRtQuantiles"))

#' @rdname chromatographyDuration
#' @export
setGeneric("chromatographyDuration", function(object, ...) standardGeneric("chromatographyDuration"))

#' @rdname rtAcquisitionRange
#' @export
setGeneric("rtAcquisitionRange", function(object, ...) standardGeneric("rtAcquisitionRange"))
