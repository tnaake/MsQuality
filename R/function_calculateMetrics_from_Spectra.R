#' @name calculateMetricsFromOneSampleSpectra
#' 
#' @title Calculate QC metrics from a Spectra object containing only spectral 
#' data from one sample
#' 
#' @description
#' The function \code{calculateMetricsFromOneSampleSpectra} calculates quality 
#' metrics from a \code{Spectra} containing spectral data from one sample. 
#' 
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments 
#' passed to the quality metric functions can be specified by the \code{params}
#' argument. \code{params} can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' The \code{Spectra} object will only contain spectral data from one 
#' data origin (\code{spectra$dataOrigin} is of length 1).
#' 
#' @param spectra \code{Spectra} object
#' @param metrics `\code{character} specifying the quality metrics to be 
#' calculated on \code{spectra}
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return named \code{numeric} vector
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @importFrom methods is
#' @importFrom Spectra Spectra
#' 
#' @examples
#' library(msdata)
#' library(Spectra)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)[1]
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#'     
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (MsLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' MsQuality:::calculateMetricsFromOneSampleSpectra(spectra = spectra, 
#'     metrics = metrics, msLevel = 1, change = "jump", relativeTo = "Q1")
#' MsQuality:::calculateMetricsFromOneSampleSpectra(spectra = spectra, 
#'     metrics = metrics, msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetricsFromOneSampleSpectra <- function(spectra, 
    metrics = qualityMetrics(spectra), ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)

    if(!is(spectra, "Spectra")) stop("'spectra' is not of class 'Spectra'")
    if(length(unique(Spectra::dataOrigin(spectra))) != 1) 
        stop("'spectra' should only contain data from one origin")
    dots <- list(...)

    ## prepare the argument for the metric functions by writing spectra to a 
    ## list
    sp_l <- list(spectra = spectra)
    args <- c(sp_l, dots)

    ## calculate the metrics (using all metrics defined in metrics) using the
    ## spectra object
    ## lapply is the outer loop that iterates through the functions `metrics`
    metrics_vals <- lapply(seq_along(metrics), function(i) {
        do.call(metrics[i], args)
    })
    
    ## add attributes 
    names(metrics_vals) <- metrics
    metrics_vals <- unlist(metrics_vals)
    attributes(metrics_vals) <- c(attributes(metrics_vals), dots)

    ## return the object
    metrics_vals
}

#' @name calculateMetricsFromSpectra
#' 
#' @title Calculate QC metrics from a Spectra object
#' 
#' @description
#' The function \code{calculateMetricsFromSpectra} calculates quality metrics 
#' from a \code{Spectra} object. The function will calculate the 
#' metrics per sample according to the \code{dataOrigin} information.
#' 
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments 
#' passed to the quality metric functions can be specified by the \code{params}
#' argument. \code{params} can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param spectra \code{Spectra} object
#' @param metrics \code{character} specifying the quality metrics to be 
#' calculated on \code{spectra}
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra (in rows)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom Spectra Spectra
#' @importFrom ProtGenerics spectra
#' @importFrom methods is
#' 
#' @examples 
#' library(msdata)
#' library(Spectra)
#' 
#' ## define file names containing spectra data for the samples
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' 
#' ## import the data and add it to the spectra object
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#' 
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' calculateMetricsFromSpectra(spectra = spectra, metrics = metrics, 
#'     msLevel = 1, change = "jump", relativeTo = "Q1")
#' calculateMetricsFromSpectra(spectra = spectra, metrics = metrics, 
#'     msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetricsFromSpectra <- function(spectra, 
    metrics = qualityMetrics(spectra), ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)
    
    if(!is(spectra, "Spectra")) stop("spectra is not of class 'Spectra'")
    
    ## get first the number of spectra in the mse object, one spectra should 
    ## refer to one mzML file/sample 
    dO <- Spectra::dataOrigin(spectra)
    dO_unique <- unique(dO)
    
    ## iterate through the different spectra per dataOrigin and calculate the 
    ## quality metrics using the calculateMetricsFromOneSampleSpectra
    ## the lapply loop returns list containing named numeric vectors
    spectra_metrics <- lapply(dO_unique, function(dO_unique_i) {
        spectra_i <- spectra[dO == dO_unique_i, ]
        calculateMetricsFromOneSampleSpectra(spectra = spectra_i, 
            metrics = metrics, ...)
    })
    df <- do.call("rbind", spectra_metrics)
    rownames(df) <- dO_unique
    
    ## add attributes
    dots <- list(...)
    attributes(df) <- c(attributes(df), dots)
    
    ## return the data.frame
    df
}

#' @name calculateMetrics
#' 
#' @title Calculate QC metrics from a Spectra object
#' 
#' @description
#' Calculate QC metrics from a `Spectra` object. 
#' `calculateMetrics` is a wrapper for the function
#' `calculateMetricsFromSpectra`.
#' 
#' @details
#' The metrics are defined by the argument `metrics`. Further arguments 
#' passed to the quality metric functions can be specified by the `params`
#' argument. `params` can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param object `Spectra` object
#' @param metrics `character` specifying the quality metrics to be calculated
#' on `object`
#' @param ... arguments passed to the quality metrics functions defined in 
#' `metrics`
#' 
#' @return `data.frame` containing in the columns the metrics for the 
#' different spectra (samples in row)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom methods is
#' @importFrom Spectra Spectra
#' 
#' @examples
#' library(msdata)
#' library(Spectra)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#'     
#' #' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (MsLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' calculateMetrics(object = spectra, metrics = metrics, 
#'     msLevel = 1, change = "jump", relativeTo = "Q1")
#' calculateMetrics(object = spectra, metrics = metrics, 
#'     msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetrics <- function(object, 
        metrics = qualityMetrics(object), ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(object), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(object), 
        several.ok = TRUE)
    
    if (is(object, "Spectra")) {
        metrics_vals <- calculateMetricsFromSpectra(spectra = object, 
            metrics = metrics, ...)
    }
    
    ## TODO: create if for MsExperiment class (to be included once MsExperiment
    ## package is accepted in BioC)
    
    ## return the object
    metrics_vals
}

