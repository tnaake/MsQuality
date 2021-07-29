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
#' @param params `list` containing parameters passed to the quality metrics
#' functions defined in `metrics`
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
    if (any(duplicated(names(params)))) stop("params contains duplicated names")
    
    ## prepare the argument for the metric functions by writing spectra to a 
    ## list
    sp_l <- list(spectra = spectra)
    
    ## calculate the metrics (using all metrics defined in metrics) using the
    ## spectra object
    ## lapply is the outer loop going iterating through the functions `metrics`
    metrics_vals <- lapply(seq_along(metrics), function(i) {
        
        ##formalArgs_i <- methods::formalArgs(metrics[i])
        formals_i <- formals(metrics[i])
        
        ## how to deal with different arguments for some of the functions? (1-4)
        ## use all parameter combinations from params that fit to the function:
        ## 1) check if parameters in params are in formals of the metric 
        ## function, update the list entry for these paramters
        params_i <- params[names(params) %in% names(formals_i)]
        inds <- match(names(formals_i), names(params_i))
        inds <- inds[!is.na(inds)]
        formals_i[names(formals_i) %in% names(params_i)] <- params_i[inds]
        ## 2) when there are calls/language types in formals_i, i.e. if there
        ## are several options for the arguments defined, take only the first 
        ## option, e.g. if we have function(a = c(1:3)) ..., we will only 
        ## continue with a = 1, NB: this is not the case if we have specified
        ## the arguments within params
        formals_i <- lapply(formals_i, function(x) 
            if (is.call(x)) {eval(x)[1]} else {x})
        ## 3) remove the type of arguments that are refObject, this will be 
        ## for instance the spectra argument
        inds_remove <- unlist(lapply(formals_i, function(x) is(x, "refObject")))
        formals_i <- formals_i[!inds_remove]
        ## 4) use all parameter combinations defined in formals_i and create a 
        ## grid
        params_i <- expand.grid(formals_i, stringsAsFactors = FALSE)
        
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
            }, simplify = FALSE)
            
            ## iterate trough the parameter combinations and return the names
            names_metric <- apply(params_i, 1, function(j) {
                paste(c(metrics[i], paste0(names(j), j)), collapse = "_")
            })
        }
        
        if (is(metric_i_j, "list")) {
            names_metric <- lapply(seq_along(names_metric), function(j) {
                if (length(names(metric_i_j[[j]])) > 0) {
                    paste(names_metric[j], names(metric_i_j[[j]]), sep = "_")
                } else {
                    paste(names_metric[j])
                }
            })
            names_metric <- unlist(names_metric)
            metric_i_j <- unlist(metric_i_j)
            
        }
        # if (is(metric_i_j, "vector") & length(metric_i_j) > 1) {
        #     names_metric <- paste(names_metric, names(metric_i_j), sep = "_")
        # }
        names(metric_i_j) <- names_metric

        return(metric_i_j)
    })
    unlist(metrics_vals)
}

#' @name calculateMetricsFromMSE
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
#' @param msexp `MsExperiment` object
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
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' calculateMetricsFromMsExperiment(msexp = mse, metrics = metrics, 
#'     params = params_l)
calculateMetricsFromMsExperiment <- function(msexp, 
    metrics = qualityMetrics(msexp), params = list()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(mse), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(mse)
    metrics <- match.arg(metrics, choices = qualityMetrics(msexp), 
                         several.ok = TRUE)
    
    if(!is(msexp, "MsExperiment")) stop("mse is not of class 'MsExperiment'")
    
    ## get first the number of spectra in the mse object, one spectra should 
    ## refer to one mzML file/sample 
    sD <- MsExperiment::sampleData(msexp)
    nsample <- nrow(sD)
    
    ## iterate through the different spectra in mse and calculate the 
    ## quality metrics using the calculateMetricsFromSpectra
    ## the lapply loop returns list containing named numeric vectors
    mse_metrics <- lapply(seq_len(nsample), function(i) {
        spectra_i <- ProtGenerics::spectra(msexp[, i])
        calculateMetricsFromSpectra(spectra = spectra_i, metrics = metrics, 
                                    params = params)
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
#' @param params `list` containing parameters passed to the quality metrics
#' functions defined in `metrics`
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
#' metrics <- c("areaUnderTIC", "rtDuration", "msSignal10XChange")
#' 
#' ## additional parameters passed to the quality metrics functions
#' ## (MSLevel is an argument of areaUnderTIC and msSignal10XChange,
#' ## relativeTo is an argument of msSignal10XChange)
#' params_l <- list(MSLevel = 1, relativeTo = c("Q1", "previous"), 
#'     change = c("jump", "fall"))
#'     
#' ## calculate the metrics
#' calculateMetrics(object = spectra, metrics = metrics, 
#'     params = params_l)
calculateMetrics <- function(object, 
        metrics = qualityMetrics(object), params = list()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(object), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(object), 
        several.ok = TRUE)
    
    if (is(object, "Spectra")) {
        metrics_vals <- calculateMetricsFromSpectra(spectra = object, 
            metrics = metrics, params = params)
    }
    
    if (is(object, "MsExperiment")) {
        metrics_vals <- calculateMetricsFromMsExperiment(msexp = object,
            metrics = metrics, params = params)
    }
    
    metrics_vals
    
}

