#' @name calculateMetricsFromOneSampleChromatograms
#'
#' @title Calculate QC metrics from a Chromatograms object containing only
#' chromatographic data from one sample
#'
#' @description
#' The function \code{calculateMetricsFromOneSampleChromatograms} calculates
#' quality metrics from a \code{Chromatograms} object containing chromatographic
#' data from one sample.
#'
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments
#' passed to the quality metric functions can be specified by the \code{params}
#' argument. \code{params} can contain named entries which are matched against
#' the formal arguments of the quality metric functions.
#'
#' The \code{Chromatograms} object will only contain chromatographic data from
#' one data origin (e.g. \code{object$dataOrigin} is of length 1). The
#' grouping is specified by the argument \code{f}.
#'
#' Setting the argument \code{filterEmptyObject} to \code{TRUE} will
#' remove zero-length entries and entries with intensities that are \code{Inf}
#' from the \code{Chromatograms} object.
#'
#' @param chromatograms \code{Chromatograms} object
#' @param metrics \code{character} specifying the quality metrics to be
#' calculated on \code{chromatograms}
#' @param filterEmptyObject \code{logical(1)} specifying if empty entries
#' and entries with intensity zero or \code{Inf} of the \code{Chromatograms}
#' object will be removed
#' @param f \code{character}, grouping parameter for \code{chromatograms}.
#'     Defaults to \code{chromatograms$dataOrigin}; if missing, \code{NULL}, or
#'     all \code{NA}, all chromatograms are treated as one sample.
#' @param ... arguments passed to the quality metrics functions defined in
#' \code{metrics}
#'
#' @return named \code{numeric} vector
#'
#' @author Philippine Louail
#'
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' library(Chromatograms)
#' ## Create a Chromatograms object with ChromBackendMemory
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem2", "mem3")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(3.5, 4.0, 4.5),
#'                intensity = c(80, 120, 90)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' metrics <- c("chromatogramDuration", "maxIntensity", "intensityMean")
calculateMetricsFromOneSampleChromatograms <- function(
    chromatograms,
    metrics = qualityMetrics(chromatograms),
    filterEmptyObject = FALSE,
    f = chromatograms$dataOrigin,
    ...
) {
    metrics <- match.arg(metrics, choices = qualityMetrics(chromatograms),
                         several.ok = TRUE)

    if (length(filterEmptyObject) != 1 | !is.logical(filterEmptyObject))
        stop("'filterEmptyObject' has to be either TRUE or FALSE")
    if (!is(chromatograms, "Chromatograms"))
        stop("'chromatograms' is not of class 'Chromatograms'")
    if (length(unique(f)) != 1)
        stop("'chromatograms' should only contain data from one origin")

    if (filterEmptyObject) {
        int <- intensity(chromatograms)
        keep <- vapply(
            int,
            function(x) {
                !is.null(x) &&
                    length(x) > 0 &&
                    !all(is.na(x)) &&
                    !all(x == 0) &&
                    !any(is.infinite(x))
            },
            logical(1)
        )
        chromatograms <- chromatograms[keep]
    }

    metrics_vals <- lapply(metrics, function(metric_name) {
        result <- get(metric_name)(chromatograms, ...)
        result_attrs <- attributes(result)
        result_names <- names(result)
        result <- unname(result)
        if (length(result) > 1 && !is.null(result_names))
            names(result) <- result_names
        for (attr_name in setdiff(names(result_attrs), "names"))
            attr(result, attr_name) <- result_attrs[[attr_name]]
        result
    })
    names(metrics_vals) <- metrics
    metrics_vals_attributes <- unlist(lapply(metrics_vals, attributes)[[1]])
    metrics_vals <- unlist(metrics_vals, use.names = TRUE)
    attributes(metrics_vals) <- c(attributes(metrics_vals), 
                                   metrics_vals_attributes, list(...))
    metrics_vals
}

#' @name calculateMetricsFromChromatograms
#'
#' @title Calculate QC metrics from a Chromatograms object
#'
#' @description
#' The function \code{calculateMetricsFromChromatograms} calculates quality
#' metrics from a \code{Chromatograms} object. The function will calculate the
#' metrics per sample according to the grouping parameter \code{f},
#' e.g. \code{dataOrigin} information.
#'
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments
#' passed to the quality metric functions can be specified by \code{...}.
#' The additional arguments \code{...} are matched against
#' the formal arguments of the quality metric functions.
#' Samples will be processed in parallel using the default parallel processing
#' setup ([bpparam()]) or with the parallel processing setup defined with
#' parameter \code{BPPARAM}.
#'
#' Setting the argument \code{filterEmptyObject} to \code{TRUE} will
#' remove zero-length entries, zero-intensity entries, and entries with
#'
#' @param chromatograms \code{Chromatograms} object
#' @param metrics \code{character} specifying the quality metrics to be
#' calculated on \code{chromatograms}
#' @param filterEmptyObject \code{logical(1)} specifying if empty entries
#' and entries with intensity zero of the \code{Chromatograms} object will be
#' removed
#' @param f \code{character} defining which chromatograms in \code{chromatograms}
#' belong to one sample. Defaults to \code{f = dataOrigin(chromatograms)}.
#' Chromatograms from the same original data file are processed together
#' (and in parallel for different files).
#' @param format \code{character(1)} output format. Only \code{"data.frame"}
#'     is currently supported.
#' @param BPPARAM Parallel processing setup. Defaults to \code{BPPARAM = bpparam()}.
#'     See [bpparam()] for details on parallel processing with \code{BiocParallel}.
#' @param ... arguments passed to the quality metrics functions defined in
#' \code{metrics}
#'
#' @return
#' A \code{data.frame} containing in the columns the metrics for the different
#' chromatograms of identical \code{dataOrigin{chromatograms}} (in rows).
#'
#' @author Philippine Louail
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom BiocParallel bplapply bpparam
#'
#' @examples
#' library(Chromatograms)
#' ## Create a Chromatograms object with ChromBackendMemory
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem2", "mem3")
#' )
#' pdata <- list(
#'     data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
#'                intensity = c(100, 250, 400, 300, 150)),
#'     data.frame(rtime = c(3.5, 4.0, 4.5),
#'                intensity = c(80, 120, 90)),
#'     data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
#'                intensity = c(80, 500, 1200, 600, 120))
#' )
#' chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
#' metrics <- c("chromatogramDuration", "maxIntensity", "intensityMean")
calculateMetricsFromChromatograms <- function(
    chromatograms,
    metrics,
    filterEmptyObject = FALSE,
    f = dataOrigin(chromatograms),
    format = c("data.frame", "mzQC"),
    ...,
    BPPARAM = bpparam()
) {
    metrics <- match.arg(metrics, choices = qualityMetrics(chromatograms),
                         several.ok = TRUE)
    if (length(filterEmptyObject) != 1 | !is.logical(filterEmptyObject))
        stop("'filterEmptyObject' has to be either TRUE or FALSE")
    format <- match.arg(format)
    if (format != "data.frame")
        stop("Only format = 'data.frame' is supported currently")

    if (!is(chromatograms, "Chromatograms"))
        stop("chromatograms is not of class 'Chromatograms'")

    if (is.null(f) || all(is.na(f))) f <- rep("sample", length(chromatograms))
    f_unique <- unique(f)

    chromatograms_metrics <- bplapply(f_unique, function(f_unique_i, ...) {
        calculateMetricsFromOneSampleChromatograms(
            chromatograms = chromatograms[f == f_unique_i], metrics = metrics,
            filterEmptyObject = filterEmptyObject, ...)
    }, ..., BPPARAM = BPPARAM)

    names(chromatograms_metrics) <- f_unique

    if (length(chromatograms_metrics) == 0) {
        return(data.frame())
    }

    if (format == "data.frame") {
        obj_attributes <- lapply(chromatograms_metrics, attributes)[[1]]
        obj <- do.call("rbind", chromatograms_metrics)
        col_names <- names(chromatograms_metrics[[1]])
        obj <- as.data.frame(obj)
        colnames(obj) <- col_names
        rownames(obj) <- f_unique
        dots <- list(...)
        attributes(obj) <- c(attributes(obj), obj_attributes, dots)
    }
    obj
}

#' @rdname calculateMetrics
#' @export
setMethod("calculateMetrics", "Chromatograms", function(object, metrics = qualityMetrics(object),
    filterEmptyObject = FALSE, ...) {
    metrics <- match.arg(metrics, choices = qualityMetrics(object),
        several.ok = TRUE)
    if (length(filterEmptyObject) != 1 | !is.logical(filterEmptyObject))
        stop("'filterEmptyObject' has to be either TRUE or FALSE")
    calculateMetricsFromChromatograms(chromatograms = object, metrics = metrics,
        filterEmptyObject = filterEmptyObject, ...)
})
