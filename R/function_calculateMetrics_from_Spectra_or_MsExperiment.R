#' @name calculateMetricFromMSE
#' 
#' @title Calculate QC metrics from a MsExperiment object
#' 
#' @description
#' The function `calculateMetricsFromMSE` calculates quality metrics from a
#' `MsExperiment`. Each spectra in the `mse` object should refer to one
#' mzML file/to one sample.
#' 
#' @details
#' The metrics are defined by the argument `metrics`. Further arguments 
#' passed to the quality metric functions can be specified by the `params`
#' argument. `params` can contain named entries which are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' @param mse `MsExperiment` object
#' @param metrics `character` specifying the quality metrics to be calculated
#' on `mse`
#' @param params `list` containing parameters passed to the quality metrics
#' functions defined in `metrics`
#' 
#' @return `data.frame` containing in the columns the metrics for the 
#' different spectra (in rows)
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom MsExperiment sampleData
#' @importFrom ProtGenerics spectra
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
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' calculateMetricsFromMsExperiment(mse = mse, metrics = metrics, 
#'     params = params_l)
calculateMetricsFromMsExperiment(mse, metrics = qualityMetrics(mse),
    params = list()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(mse), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(mse)
    metrics <- match.arg(metrics, choices = qualityMetrics(mse), 
        several.ok = TRUE)
    
    if(!is(mse, "MsExperiment")) stop("mse is not of class 'MsExperiment'")
    
    ## get first the number of spectra in the mse object, one spectra should 
    ## refer to one mzML file/sample 
    sD <- MsExperiment::sampleData(mse)
    nsample <- nrow(sD)
    
    ## iterate through the different spectra in mse and calculate the 
    ## quality metrics using the calculateMetricsFromSpectra
    ## the lapply loop returns list containing named numeric vectors
    mse_metrics <- lapply(seq_len(nsample), function(i) {
        spectra_i <- ProtGenerics::spectra(mse[, i])
        calculateMetricsFromSpectra(spectra = spectra_i, metrics = metrics, 
            params = params)
    })
    df <- do.call("rbind", mse_metrics)
    rownames(df) <- rownames(sD)
    return(df)
}

#' @name calculateMetricsFromSpectra
#' 
#' @title Calculate QC metrics from a Spectra object
#' 
#' @description
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
#' @param params `list` containing parameters passed to the quality metrics
#' functions defined in `metrics`
#' 
#' @return named `numeric` vector
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom methods formalArgs
#' 
#' @examples
#' library(msdata)
#' library(Spectra)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' calculateMetricsFromSpectra(spectra = spectra, metrics = metrics, 
#'     params = params_l)
calculateMetricsFromSpectra <- function(spectra, 
    metrics = qualityMetrics(spectra), params = list()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)
    
    if(!is(spectra, "Spectra")) stop("spectra is not of class 'Spectra'")
    
    ## prepare the argument for the metric functions by writing spectra to a 
    ## list
    sp_l <- list(spectra = spectra)
    
    ## calculate the metrics (using all metrics defined in metrics) using the
    ## spectra object
    ## lapply is the outer loop going iterating through the functions `metrics`
    metrics_vals <- lapply(seq_along(metrics), function(i) {
        formalArgs_i <- methods::formalArgs(metrics[i])
        ## how to deal with different arguments for some of the functions?
        ## use all parameter combinations from params that fit to the function
        params_i <- params[names(params) %in% formalArgs_i]
        params_i <- expand.grid(params_i, stringsAsFactors = FALSE)
        
        ## in case there are no further parameters for the function metrics[i]:
        if (nrow(params_i) == 0) {
            metric_i_j <- do.call(metrics[i], args = sp_l)
            names_metric <- metrics[i]
            
            ## in case there are further parameters for the function metrics[i]:
        } else {
            ## iterate through the parameter combinations and return the values
            metric_i_j <- apply(params_i, 1, function(j) {
                args_l <- append(sp_l, as.list(j))
                do.call(metrics[i], args = args_l)
            })
            
            ## iterate trough the parameter combinations and return the names
            names_metric <- apply(params_i, 1, function(j) {
                paste(c(metrics[i], paste0(names(j), j)), collapse = "_")
            })
        }
        
        if (is(metric_i_j, "matrix")) {
            names_metric <- paste(names_metric, rownames(metric_i_j), sep = "_")
            metric_i_j <- as.vector(metric_i_j)
            
        } 
        names(metric_i_j) <- names_metric
        
        
        return(metric_i_j)
    })
    metrics_vals <- unlist(metrics_vals)
    
    return(metrics_vals)
}
