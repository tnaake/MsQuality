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
#' data origin (e.g. \code{spectra$dataOrigin} is of length 1). The grouping 
#' is specified by the argument \code{f}.
#' 
#' @param spectra \code{Spectra} object
#' @param metrics \code{character} specifying the quality metrics to be 
#' calculated on \code{spectra}
#' @param f \code{character}, grouping parameter for \code{spectra}
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
    metrics = qualityMetrics(spectra), f = spectra$dataOrigin, ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)

    if(!is(spectra, "Spectra")) stop("'spectra' is not of class 'Spectra'")
    if(length(unique(f)) != 1) 
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
#' metrics per sample according to the grouping parameter \code{f}, 
#' e.g. \code{dataOrigin} information. Samples will be processed in parallel
#' using the default parallel processing setup ([bpparam()]) or with the
#' parallel processing setup defined with parameter `BPPARAM`.
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
#' @param f \code{character} defining which spectra in `spectra` belong to
#'     one sample. Defaults to `f = dataOrigin(spectra)`. Spectra from the
#'     same original data file are processed together (and in parallel for
#'     different files).
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' @param BPPARAM Parallel processing setup. Defaults to `BPPARAM = bpparam()`.
#'     See [bpparam()] for details on parallel processing with `BiocParallel`.
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra (in rows)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}, Johannes Rainer
#' 
#' @export
#' 
#' @importFrom Spectra Spectra
#' @importFrom methods is
#' @importMethodsFrom Spectra dataOrigin
#' @importFrom BiocParallel bplapply bpparam
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
    metrics = qualityMetrics(spectra), f = dataOrigin(spectra), ...,
    BPPARAM = bpparam()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)
    
    if(!is(spectra, "Spectra")) stop("spectra is not of class 'Spectra'")
    
    ## get first the number of spectra in the mse object, one spectra should 
    ## refer to one mzML file/sample or other grouping factor specified by
    ## paramter f
    f_unique <- unique(f)
    
    ## iterate through the different spectra per dataOrigin and calculate the 
    ## quality metrics using the calculateMetricsFromOneSampleSpectra
    ## the lapply loop returns list containing named numeric vectors
    spectra_metrics <- bplapply(f_unique, function(f_unique_i, ...) {
        calculateMetricsFromOneSampleSpectra(
            spectra = spectra[f == f_unique_i], metrics = metrics, ...)
    }, ..., BPPARAM = BPPARAM)
    df <- do.call("rbind", spectra_metrics)
    rownames(df) <- f_unique
    
    ## add attributes
    dots <- list(...)
    attributes(df) <- c(attributes(df), dots)
    
    ## return the data.frame
    df
}


#' @name calculateMetricsFromMsExperiment
#' 
#' @title Calculate QC metrics from a MsExperiment object
#' 
#' @description
#' The function \code{calculateMetricsFromMsExperiment} calculates quality 
#' metrics from a \code{MsExperiment} object. Each spectra in the 
#' \code{msexp} object should refer to one mzML file/to one sample.
#' 
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments 
#' passed to the quality metric functions can be specified by the \code{params}
#' argument. \code{params} can contain named entries which are matched against
#' the formal arguments of the quality metric functions.
#' 
#' @param msexp \code{MsExperiment} object
#' @param metrics \code{character} specifying the quality metrics to be 
#' calculated on \code{msexp}
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra (in rows)
#'
#' @inheritParams calculateMetricsFromSpectra
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom Spectra Spectra
#' @importFrom ProtGenerics spectra
#' @importFrom MsExperiment MsExperiment sampleData
#' @importFrom methods is
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' 
#' msexp <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(msexp) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(msexp) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' msexp <- linkSampleData(msexp, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' msexp <- linkSampleData(msexp, with = "experimentFiles.annotations",
#'      sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(msexp) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' calculateMetricsFromMsExperiment(msexp = msexp, metrics = metrics,
#'     msLevel = 1, change = "jump", relativeTo = "Q1")
#'     
#' calculateMetricsFromMsExperiment(msexp = msexp, metrics = metrics, 
#'     msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetricsFromMsExperiment <- function(msexp, 
    metrics = qualityMetrics(msexp), ..., BPPARAM = bpparam()) {
  
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(mse), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(mse)
    metrics <- match.arg(metrics, choices = qualityMetrics(msexp), 
        several.ok = TRUE)
    
    if(!is(msexp, "MsExperiment")) 
        stop("'msexp' is not of class 'MsExperiment'")
    
    
    ## get Spectra object from MsExperiment object and calculate the quality  
    ## metrics using the calculateMetricsFromSpectra function, the metrics
    ## will be stored in the data.frame df
    sps <- spectra(msexp)
    df <- calculateMetricsFromSpectra(spectra = sps, metrics = metrics, ...,
                                      BPPARAM = BPPARAM)
    
    ## return the data.frame
    df
}

#' @name calculateMetrics
#' 
#' @title Calculate QC metrics from a Spectra or MsExperiment object
#' 
#' @description
#' Calculate QC metrics from a \code{Spectra} or \code{MsExperiment} object. 
#' \code{calculateMetrics} is a wrapper for the functions
#' \code{calculateMetricsFromSpectra} and 
#' \code{calculateMetricsFromMsExperiment}.
#' 
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments 
#' passed to the quality metric functions can be specified by the \code{params}
#' argument. \code{params} can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param object \code{Spectra} or \code{MsExperiment} object
#' @param metrics \code{character} specifying the quality metrics to be 
#' calculated on \code{object}
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra and in rows the samples
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
    
    if (is(object, "MsExperiment")) {
      metrics_vals <- calculateMetricsFromMsExperiment(msexp = object, 
            metrics = metrics, ...)
    }
    
    ## return the object
    metrics_vals
}

