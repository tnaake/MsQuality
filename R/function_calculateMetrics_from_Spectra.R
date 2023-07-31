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
#' @param filterEmptySpectra \code{logical(1)} specifying if empty entries and
#' entries with intensity zero of the \code{Spectra} object will be removed
#' @param f \code{character}, grouping parameter for \code{spectra}
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return named \code{numeric} vector
#' 
#' @author Thomas Naake
#' 
#' @importFrom methods is
#' @importFrom Spectra Spectra filterIntensity filterEmptySpectra
#' @examples
#' library(msdata)
#' library(Spectra)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)[1]
#' spectra <- Spectra(fls, backend = MsBackendMzR())
#' 
#' ## define the quality metrics to be calculated
#' metrics <- c("areaUnderTic", "chromatographyDuration", "msSignal10xChange")
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
    metrics = qualityMetrics(spectra), filterEmptySpectra = FALSE, 
    f = spectra$dataOrigin, ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)
    
    if (length(filterEmptySpectra) != 1 & !is.logical(filterEmptySpectra))
        stop("'filterEmptySpectra' has to be either TRUE or FALSE")

    if(!is(spectra, "Spectra")) stop("'spectra' is not of class 'Spectra'")
    
    if(length(unique(f)) != 1) 
        stop("'spectra' should only contain data from one origin")
    
    ## in case of filterEmptySpectra == TRUE, remove the entries with
    ## zero or Inf intensity and remove the entries with empty spectra
    if (filterEmptySpectra) {
        spectra <- spectra |>
            filterIntensity(intensity = c(0, Inf)) |>
            filterEmptySpectra(spectra)
    }
    
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
    
    ## add attributes (attributes of metrics_vals and dots)
    names(metrics_vals) <- metrics
    metrics_vals_attributes <- lapply(metrics_vals, attributes)
    names(metrics_vals_attributes) <- NULL
    metrics_vals_attributes <- unlist(metrics_vals_attributes)
    
    metrics_vals <- unlist(metrics_vals)
    attributes(metrics_vals) <- c(attributes(metrics_vals), 
        metrics_vals_attributes, dots)

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
#' e.g. \code{dataOrigin} information. 
#' 
#' Two format options are available:
#' 
#' - \code{format = "data.frame"} returns the metrics as a \code{data.frame}, \cr
#' - \code{format = "mzQC"} returns the metrics as a list of \code{MzQCmzQC} 
#'   objects. \cr
#' 
#' @details
#' The metrics are defined by the argument \code{metrics}. Further arguments 
#' passed to the quality metric functions can be specified by \code{...}. 
#' The additional arguments \code{...} are matched against 
#' the formal arguments of the quality metric functions. 
#' 
#' Samples will be processed in parallel
#' using the default parallel processing setup ([bpparam()]) or with the
#' parallel processing setup defined with parameter \code{BPPARAM}.
#' 
#' @param spectra \code{Spectra} object
#' @param metrics \code{character} specifying the quality metrics to be 
#' calculated on \code{spectra}
#' @param filterEmptySpectra \code{logical(1)} specifying if empty entries and
#' entries with intensity zero of the \code{Spectra} object will be removed
#' @param f \code{character} defining which spectra in \code{spectra} belong to
#'     one sample. Defaults to \code{f = dataOrigin(spectra)}. Spectra from the
#'     same original data file are processed together (and in parallel for
#'     different files).
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
#' the columns the metrics for the different spectra of identical 
#' \code{dataOrigin{spectra}} (in rows).
#' In case of \code{format = "mzQC"}, a \code{list} of \code{MzQCmzQC} objects
#' containing the metrics for the different spectra of identical 
#' \code{dataOrigin{spectra}}
#' 
#' @author Thomas Naake, Johannes Rainer
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
#' metrics <- c("areaUnderTic", "chromatographyDuration", "msSignal10xChange")
#' 
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (msLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' 
#' ## format = "data.frame"
#' calculateMetricsFromSpectra(spectra = spectra, metrics = metrics,
#'     format = "data.frame", msLevel = 1, change = "jump", relativeTo = "Q1")
#' calculateMetricsFromSpectra(spectra = spectra, metrics = metrics, 
#'     format = "data.frame", msLevel = 1, change = "fall", 
#'     relativeTo = "previous")
#'     
#' ## format = "mzQC"
#' ##calculateMetricsFromSpectra(spectra = spectra, metrics = metrics,
#' ##    format = "mzQC", msLevel = 1, change = "jump", relativeTo = "Q1")
#' ##calculateMetricsFromSpectra(spectra = spectra, metrics = metrics, 
#' ##    format = "mzQC", msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetricsFromSpectra <- function(spectra, metrics, 
    filterEmptySpectra = FALSE, f = dataOrigin(spectra), 
    format = c("data.frame", "mzQC"), ..., BPPARAM = bpparam()) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(spectra), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(spectra), 
        several.ok = TRUE)
    
    if (length(filterEmptySpectra) != 1 & !is.logical(filterEmptySpectra))
        stop("'filterEmptySpectra' has to be either TRUE or FALSE")
    
    format <- match.arg(format)
    
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
            spectra = spectra[f == f_unique_i], metrics = metrics, 
            filterEmptySpectra = FALSE, ...)
    }, ..., BPPARAM = BPPARAM)
    
    ## add file names as names of the list
    names(spectra_metrics) <- f_unique
    
    ## if format == "data.frame"
    if (format == "data.frame") {
        obj_attributes <- lapply(spectra_metrics, attributes)[[1]]
        obj <- do.call("rbind", spectra_metrics)
        
        ## add attributes
        dots <- list(...)
        attributes(obj) <- c(attributes(obj), obj_attributes, dots)
        names(obj) <- NULL
    }
    
    if (format == "mzQC") {
        obj <- transformIntoMzQC(spectra_metrics)
    }
    
    ## return the data.frame or the list of mzQC
    obj
}

#' @name transformIntoMzQC
#' 
#' @title Transform the metrics into a list of \code{MzQCmzQC} objects
#' 
#' @description
#' The function \code{transformIntoMzQC} transfers the metrics stored in 
#' \code{spectra_metrics} into a list of \code{MzQCmzQC} objects. Each list 
#' entry will refer to the corresponding entry in \code{spectra_metrics}.
#' As such, each entry contains information from a single \code{dataOrigin}
#' of a \code{Spectra} object.
#' 
#' The function \code{transformIntoMzQC} is a helper function within
#' \code{calculateMetricsFromSpectra}.
#' 
#' @details
#' The \code{MzQCmzQC} object will only contain those quality metrics
#' that have a corresponding attribute with a [PSI:MS] identifier. The 
#' matching is done via the names of each vector in \code{spectra_metrics}.
#' 
#' The Field \code{"version"} is set to the current version of the \code{rmzqc}
#' package.
#' 
#' The entry of \code{"MzQCanalysisSoftware"} is filled with the [PSI:MS] id
#' of \code{MsQuality} ("MS:") and the version is taken from
#' \code{packageDescription("MsQuality")[["Version"]]}.
#' 
#' @param spectra_metrics list of named vector
#' 
#' @return \code{list} containing as entries \code{MzQCmzQC} objects for each
#' \code{Spectra} with same \code{dataOrigin}
#' 
#' @author Thomas Naake, Johannes Rainer
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
#' metrics <- c("areaUnderTic", "chromatographyDuration", "msSignal10xChange")
#' 
#' ## obtain the spectra_metrics object
#' f <- dataOrigin(spectra)
#' f_unique <- unique(f)
#' ## spectra_metrics <- bplapply(f_unique, function(f_unique_i) {
#' ##calculateMetricsFromOneSampleSpectra(
#' ##    spectra = spectra[f == f_unique_i], metrics = metrics)
#' ##    }, BPPARAM = bpparam())
#' 
#' ## transform into mzQC objects
#' ##transformIntoMzQC(spectra_metrics)
#' 
#' @importFrom rmzqc getCVTemplate filenameToCV toAnalysisSoftware toQCMetric
#' @importFrom rmzqc getDefaultCV
#' @importClassesFrom rmzqc MzQCrunQuality 
#' @importClassesFrom rmzqc MzQCmetadata 
#' @importClassesFrom rmzqc MzQCinputFile 
#' @importClassesFrom rmzqc MzQCmzQC
#' @importClassesFrom rmzqc MzQCDateTime
#' @importFrom utils packageDescription
transformIntoMzQC <- function(spectra_metrics) {
    
    ## create mzQC objects per sample and return as a list
    res <- lapply(seq_along(spectra_metrics), function(i) {
        
        ## obtain raw file and file format
        raw_file <- names(spectra_metrics)[i]
        file_format <- getCVTemplate(accession = filenameToCV(raw_file))
        
        ## obtain information on the MsQuality package
        software <- toAnalysisSoftware(id = "MS:4000151", 
            version = packageDescription("MsQuality")$Version)
        
        ## obtain information on the run qualities
        ## find first all metrics with "MS:NNNNNNN" attributes
        spectra_metrics_i <- spectra_metrics[[i]]
        spectra_metrics_names_i <- lapply(
            strsplit(names(spectra_metrics_i), split = "[.]"), "[", 1) |>
            unlist()
        names(spectra_metrics_i) <- spectra_metrics_names_i
        attributes_i <- attributes(spectra_metrics_i)[
            names(attributes(spectra_metrics_i)) %in% spectra_metrics_names_i]
        
        ## iterate through all valid attributes, obtain the value of the metric,
        ## and rename the entry
        qc_metric_i <- lapply(seq_along(attributes_i)[1], function(j) { #######
            id_j <- attr(x = spectra_metrics_i, 
                which = names(attributes_i)[j], exact = TRUE)
            value_j <- spectra_metrics_i[names(attributes_i)[j]] |>
                as.numeric()
            toQCMetric(id = id_j, value = value_j)
        })
        
        
        run_qc <- MzQCrunQuality(
            metadata = MzQCmetadata(
                label = raw_file,
                inputFiles = list(MzQCinputFile(
                    basename(raw_file), raw_file, file_format)),
                analysisSoftware = list(software)),
            qualityMetrics = qc_metric_i
        )
        
        ## create the final object and return
        MzQCmzQC(
            version = packageDescription("rmzqc")$Version,
            creationDate = MzQCDateTime(), 
            contactName = Sys.info()[["user"]], 
            #contactAddress = "test@user.info", 
            description = paste("A mzQC document on the sample", basename(raw_file)),
            runQualities = list(run_qc),
            setQualities = list(), 
            controlledVocabularies = list(getDefaultCV()))
        
        
    })
    
    res
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
#' @param filterEmptySpectra \code{logical(1)} specifying if empty entries and
#' entries with intensity zero of the \code{Spectra} object will be removed
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra (in rows)
#'
#' @inheritParams calculateMetricsFromSpectra
#' 
#' @author Thomas Naake
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
#' metrics <- c("areaUnderTic", "chromatographyDuration", "msSignal10xChange")
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
    metrics = qualityMetrics(msexp), filterEmptySpectra = FALSE,
    ..., BPPARAM = bpparam()) {
  
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(mse), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(mse)
    metrics <- match.arg(metrics, choices = qualityMetrics(msexp), 
        several.ok = TRUE)
    
    if (length(filterEmptySpectra) != 1 & !is.logical(filterEmptySpectra))
        stop("'filterEmptySpectra' has to be either TRUE or FALSE")
    
    if(!is(msexp, "MsExperiment")) 
        stop("'msexp' is not of class 'MsExperiment'")
    
    ## get Spectra object from MsExperiment object and calculate the quality  
    ## metrics using the calculateMetricsFromSpectra function, the metrics
    ## will be stored in the data.frame df
    sps <- spectra(msexp)
    res <- calculateMetricsFromSpectra(spectra = sps, metrics = metrics, 
        filterEmptySpectra = filterEmptySpectra, ...,
        BPPARAM = BPPARAM)
    
    ## return the object
    res
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
#' @param filterEmptySpectra \code{logical(1)} specifying if empty entries and
#' entries with intensity zero of the \code{Spectra} object will be removed
#' @param ... arguments passed to the quality metrics functions defined in 
#' \code{metrics}
#' 
#' @return \code{data.frame} containing in the columns the metrics for the 
#' different spectra and in rows the samples
#' 
#' @author Thomas Naake
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
#' metrics <- c("areaUnderTic", "chromatographyDuration", "msSignal10xChange")
#'     
#' ## calculate the metrics
#' ## additional parameters passed to the quality metrics functions
#' ## (MsLevel is an argument of areaUnderTic and msSignal10xChange,
#' ## relativeTo is an argument of msSignal10xChange) passed to ...
#' calculateMetrics(object = spectra, metrics = metrics, 
#'     msLevel = 1, change = "jump", relativeTo = "Q1")
#' calculateMetrics(object = spectra, metrics = metrics, 
#'     msLevel = 1, change = "fall", relativeTo = "previous")
calculateMetrics <- function(object, 
        metrics = qualityMetrics(object), filterEmptySpectra = FALSE, 
        ...) {
    
    ## match metrics against the possible quality metrics defined in 
    ## qualityMetrics(object), throw an error if there are metrics that 
    ## are not defined in qualityMetrics(spectra)
    metrics <- match.arg(metrics, choices = qualityMetrics(object), 
        several.ok = TRUE)
    
    if (length(filterEmptySpectra) != 1 & !is.logical(filterEmptySpectra))
        stop("'filterEmptySpectra' has to be either TRUE or FALSE")
    
    if (is(object, "Spectra")) {
        metrics_vals <- calculateMetricsFromSpectra(spectra = object, 
            metrics = metrics, filterEmptySpectra = filterEmptySpectra, ...)
    }
    
    if (is(object, "MsExperiment")) {
      metrics_vals <- calculateMetricsFromMsExperiment(msexp = object, 
            metrics = metrics, filterEmptySpectra = filterEmptySpectra, ...)
    }
    
    ## return the object
    metrics_vals
}

