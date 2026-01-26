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
