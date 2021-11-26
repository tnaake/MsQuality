#' @name calculateMetricsFromSpectra
#' 
#' @title Calculate QC metrics from a Spectra object
#' 
#' @description
#' The function `calculateMetricsFromSpectra` calculates quality metrics from a
#' `Spectra`. 
#' 
#' @details
#' The metrics are defined by the argument `metrics`. Further arguments 
#' passed to the quality metric functions can be specified by the `params`
#' argument. `params` can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param spectra `Spectra` object
#' @param metrics `character` specifying the quality metrics to be calculated
#' on `mse`
#' @param ... arguments passed to the quality metrics functions defined in 
#' `metrics`
#' 
#' @return named `numeric` vector
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom methods is
#' @importFrom Spectra Spectra
#' @importFrom MsExperiment MsExperiment
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
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (MsLevel is an argument of areaUnderTic and msSignal10xChange,
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
    ##if (any(duplicated(names(params)))) stop("params contains duplicated names")
    
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
    names(metrics_vals) <- metrics
    metrics_vals <- unlist(metrics_vals)
    attributes(metrics_vals) <- c(attributes(metrics_vals), dots)
    metrics_vals
    
}

#' @name calculateMetricsFromMsExperiment
#' 
#' @title Calculate QC metrics from a MsExperiment object
#' 
#' @description
#' The function `calculateMetricsFromMsExperiment` calculates quality metrics 
#' from a `MsExperiment` object. Each spectra in the `mse` object should 
#' refer to one mzML file/to one sample.
#' 
#' @details
#' The metrics are defined by the argument `metrics`. Further arguments 
#' passed to the quality metric functions can be specified by the `params`
#' argument. `params` can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param msexp `MsExperiment` object
#' @param metrics `character` specifying the quality metrics to be calculated
#' on `mse`
#' @param ... arguments passed to the quality metrics functions defined in 
#' `metrics`
#' 
#' @return `data.frame` containing in the columns the metrics for the 
#' different spectra (in rows)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom Spectra Spectra
#' @importFrom MsExperiment MsExperiment sampleData
#' @importFrom ProtGenerics spectra
#' @importFrom methods is
#' 
#' @examples 
#' library(msdata)
#' library(MsExperiment)
#' library(S4Vectors)
#' mse <- MsExperiment()
#' sd <- DataFrame(sample_id = c("QC1", "QC2"),
#'     sample_name = c("QC Pool", "QC Pool"), injection_idx = c(1, 3))
#' sampleData(mse) <- sd
#' 
#' ## define file names containing spectra data for the samples and
#' ## add them, along with other arbitrary files to the experiment
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' experimentFiles(mse) <- MsExperimentFiles(
#'     mzML_files = fls,
#'     annotations = "internal_standards.txt")
#' ## link samples to data files: first sample to first file in "mzML_files",
#' ## second sample to second file in "mzML_files"
#' mse <- linkSampleData(mse, with = "experimentFiles.mzML_files",
#'     sampleIndex = c(1, 2), withIndex = c(1, 2))
#' mse <- linkSampleData(mse, with = "experimentFiles.annotations",
#'                       sampleIndex = c(1, 2), withIndex = c(1, 1))
#'
#' library(Spectra)
#' ## import the data and add it to the mse object
#' spectra(mse) <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "rtDuration", "msSignal10xChange")
#' 
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     msLevel = 1, change = "jump", relativeTo = "Q1")
#' calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetricsFromMsExperiment <- function(msexp, 
    metrics = qualityMetrics(msexp), ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(mse), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(mse)
    metrics <- match.arg(metrics, choices = qualityMetrics(msexp), 
                         several.ok = TRUE)
    
    if(!is(msexp, "MsExperiment")) stop("mse is not of class 'MsExperiment'")
    
    ## get first the number of spectra in the mse object, one spectra should 
    ## refer to one mzML file/sample 
    sD <- sampleData(msexp)
    nsample <- nrow(sD)
    
    ## iterate through the different spectra in mse and calculate the 
    ## quality metrics using the calculateMetricsFromSpectra
    ## the lapply loop returns list containing named numeric vectors
    mse_metrics <- lapply(seq_len(nsample), function(i) {
        spectra_i <- spectra(msexp[, i])
        calculateMetricsFromSpectra(spectra = spectra_i, 
            metrics = metrics, ...)
    })
    df <- do.call("rbind", mse_metrics)
    rownames(df) <- rownames(sD)
    df
}

#' @name calculateMetrics
#' 
#' @title Calculate QC metrics from a Spectra or MsExperiment object
#' 
#' @description
#' Calculate QC metrics from a `Spectra` or `MsExperiment` object. 
#' `calculateMetrics` is a wrapper for the functions
#' `calculateMetricsFromSpectra` and `calculateMetricsFromMSE`
#' 
#' @details
#' The metrics are defined by the argument `metrics`. Further arguments 
#' passed to the quality metric functions can be specified by the `params`
#' argument. `params` can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param object `Spectra` or `MsExperiment` object
#' @param metrics `character` specifying the quality metrics to be calculated
#' on `mse`
#' @param ... arguments passed to the quality metrics functions defined in 
#' `metrics`
#' 
#' @return named `numeric` vector (if `object` is a `Spectra` object) or 
#' `data.frame` containing in the columns the metrics for the 
#' different spectra (in rows, if `object` is a `MsExperiment` object)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom methods is
#' @importFrom Spectra Spectra
#' @importFrom MsExperiment MsExperiment
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
    
    metrics_vals
    
}

