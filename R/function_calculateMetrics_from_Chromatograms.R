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
#' Each metric returns a per-chromatogram result: scalar metrics return a
#' \code{numeric} vector of \code{length(chromatograms)}, while multi-value
#' metrics return a \code{matrix} with \code{nrow = length(chromatograms)}.
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
#' @return named \code{list}. Each element is named by the metric and contains
#' either a \code{numeric} vector or a \code{matrix} (for multi-value metrics).
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
        get(metric_name)(chromatograms, ...)
    })
    names(metrics_vals) <- metrics
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
#' intensities that are \code{Inf} from the \code{Chromatograms} object.
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
#' @param format \code{character(1)} specifying if metrics are returned
#' as a \code{data.frame} (\code{format = "data.frame"}) or as a list of
#' \code{MzQCmzQC} objects (\code{format = "mzQC"})
#' @param BPPARAM Parallel processing setup. Defaults to \code{BPPARAM = bpparam()}.
#'     See [bpparam()] for details on parallel processing with \code{BiocParallel}.
#' @param ... arguments passed to the quality metrics functions defined in
#' \code{metrics}
#'
#' @return
#' In case of \code{format = "data.frame"}, a \code{data.frame} containing in
#' the columns the metrics for the different chromatograms of identical
#' \code{dataOrigin{chromatograms}} (in rows).
#' In case of \code{format = "mzQC"}, a \code{list} of \code{MzQCmzQC} objects
#' containing the metrics for the different chromatograms of identical
#' \code{dataOrigin{chromatograms}}
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
        ## Convert each sample's named list into a data.frame row-block.
        ## Scalar metrics become columns; matrix metrics get one column per
        ## original column, named "metric.colname".
        sample_dfs <- lapply(chromatograms_metrics, function(sample_res) {
            cols <- list()
            for (nm in names(sample_res)) {
                val <- sample_res[[nm]]
                if (is.matrix(val)) {
                    for (cn in colnames(val)) {
                        cols[[paste0(nm, ".", cn)]] <- val[, cn]
                    }
                } else {
                    cols[[nm]] <- val
                }
            }
            ## All vectors should be the same length (= number of chroms
            ## in this sample).
            as.data.frame(cols, check.names = FALSE)
        })
        obj <- do.call(rbind, sample_dfs)
    }

    if (format == "mzQC") {
        ## Convert per-sample Chromatograms metrics (named lists of
        ## vectors/matrices) into the named-numeric-vector format that
        ## transformIntoMzQC expects.
        chromatograms_metrics_flat <- lapply(
            chromatograms_metrics, function(sample_res) {
                vals <- c()
                all_attrs <- list()
                for (nm in names(sample_res)) {
                    val <- sample_res[[nm]]
                    metric_attr <- attr(val, nm)
                    if (is.matrix(val)) {
                        for (cn in colnames(val)) {
                            col_nm <- paste0(nm, ".", cn)
                            col_vals <- val[, cn]
                            names(col_vals) <- rep(col_nm, length(col_vals))
                            vals <- c(vals, col_vals)
                        }
                    } else {
                        names(val) <- rep(nm, length(val))
                        vals <- c(vals, val)
                    }
                    if (!is.null(metric_attr))
                        all_attrs[[nm]] <- metric_attr
                }
                attributes(vals) <- c(attributes(vals), all_attrs)
                vals
            })
        names(chromatograms_metrics_flat) <- names(chromatograms_metrics)
        obj <- transformIntoMzQC(chromatograms_metrics_flat)
    }

    obj
}


