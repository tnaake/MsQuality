#' @name rtDuration
#' 
#' @title RT duration (QC:4000053)
#' 
#' @description
#' "The retention time duration of the MS run in seconds, similar to the 
#' highest scan time minus the lowest scan time." [PSI:QC]
#' id: QC:4000053
#' 
#' The metric is calculated as follows:
#' (1) the retention time associated to the individual Spectra is obtained,
#' 
#' (2) the maximum and the minimum of the retention time is obtained,
#' 
#' (3) the difference between the maximum and the minimum is calculated and 
#' returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' 
#' Retention time values that are `NA` are removed.
#' 
#' @param spectra `Spectra` object
#' 
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importMethodsFrom Spectra rtime
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' rtDuration(spectra = sps)
rtDuration <- function(spectra, ...) {  
    RT <- rtime(object = spectra)
    max(RT, na.rm = TRUE) - min(RT, na.rm = TRUE)
}

#' @name rtOverTicQuantile
#' 
#' @title RT over TIC quantile (QC:4000054)
#' 
#' @description 
#' "The interval when the respective quantile of the TIC accumulates divided by 
#' retention time duration. The number of quantiles observed is given by the 
#' size of the tuple." [PSI:QC]
#' id: QC:4000054
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is ordered according to the retention time,
#' 
#' (2) the cumulative sum of the ion count is calculated (TIC),
#' 
#' (3) the quantiles are calculated according to the `probs` argument, e.g.
#' when `probs` is set to `c(0, 0.25, 0.5, 0.75, 1)` the 0%, 25%, 50%, 75% and
#' 100% quantile is calculated,
#' 
#' (4) the retention time/relative retention time (retention time divided by 
#' the total run time taking into account the minimum retention time) is 
#' calculated,
#' 
#' (5) the (relative) duration of the LC run after which the cumulative
## TIC exceeds (for the first time) the respective quantile of the
## cumulative TIC is calculated and returned.
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' @param spectra `Spectra` object
#'
#' @param probs `numeric` defining the quantiles. See [stats::quantile()] for
#'     details. Defaults to `probs = seq(0, 1, 0.25)`.
#' 
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @param relative `logical`, if set to `TRUE` the relative retention time 
#' will be returned instead of the abolute retention time
#'
#' @param ... additional arguments passed to [stats::quantile()].
#' 
#' @return `numeric` of length equal to length `probs` with the relative
#'    duration (duration divided by the total run time) after which the TIC
#'    exceeds the respective quantile of the TIC.
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}, Johannes Rainer
#' 
#' @export 
#' 
#' @importMethodsFrom Spectra ionCount filterMsLevel
#'
#' @importFrom stats quantile
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' rtOverTicQuantile(spectra = sps, msLevel = 2L)
rtOverTicQuantile <- function(spectra, probs = seq(0, 1, 0.25),# na.rm = FALSE,
    msLevel = 1L, relative = TRUE, ...) {
    
    ## truncate spectra based on the MS level
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    RT <- rtime(spectra)

    ## what is the (relative) duration of the LC run after which the cumulative
    ## TIC exceeds (for the first time) the respective quantile of the
    ## cumulative TIC? The "accumulates" is interpreted as the 
    ## "sum of the TICs of all previous spectra".
    TIC <- cumsum(ionCount(spectra))
    ticQuantile <- quantile(TIC, probs = probs, ...)
    idxs <- unlist(lapply(ticQuantile, function(z) which.max(TIC >= z)))
    
    if (relative) {
        rtMin <- min(RT)
        rtDuration <- rtDuration(spectra)
        res <- (RT[idxs] - rtMin) / rtDuration  
    } else {
        res <- RT[idxs]
    }
    
    names(res) <- names(idxs)
    res
}

#' @title Order Spectra according to increasing retention time
#' 
#' @description 
#' The function `.rtOrderSpectra` orders the features in a `Spectra` object
#' according to the (increasing) retention time values. 
#' 
#' @details
#' Internal function in quality metric functions.
#' 
#' @param spectra `Spectra` object
#' 
#' @return `Spectra` object with the features ordered according to the 
#' (increasing) retention time
#' 
#' @author Johannes Rainer
#' 
#' @importFrom ProtGenerics rtime
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0001847", "HMDB0000001", "HMDB0000001"),
#'     name = c("Caffeine", "1-Methylhistidine", "1-Methylhistidine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
#' spd$intensity <- list(
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
#' spd$rtime <- c(15.84, 9.44, 9.44)
#' sps <- Spectra(spd)
#' .rtOrderSpectra(sps)
.rtOrderSpectra <- function(spectra) {
    RT <- rtime(spectra)
    if (any(is.na(RT)))
        warning("missing retention time values. Will keep original ",
                "ordering of spectra.")
    else if (is.unsorted(RT))
        spectra <- spectra[order(RT)]
    spectra
}

#' @name rtOverMsQuarters
#' 
#' @title MS1/MS2 quantiles RT fraction (QC:4000055/QC:4000056)
#' 
#' @description
#' "The interval used for acquisition of the first, second, third, and fourth 
#' quarter of all MS1 events divided by RT-Duration." [PSI:QC]
#' id: QC:4000055
#'
#' "The interval used for acquisition of the first, second, third, and fourth 
#' quarter of all MS2 events divided by RT-Duration." [PSI:QC]
#' id: QC:4000056
#' 
#' The metric is calculated as follows:
#' (1) the retention time duration of the whole `Spectra` object is determined
#' (taking into account all the MS levels),
#' 
#' (1) the `Spectra` object is filtered according to the MS level and 
#' subsequently ordered according to the retention time
#' 
#' (2) the MS events are split into four (approximately) equal parts,
#' 
#' (3) the relative retention time is calculated (using the retention time 
#' duration from (1) and taking into account the minimum retention time),
#' 
#' (4) the relative retention time values associated to the MS event parts
#' are returned.
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' is_a: QC:4000023 ! MS1 metric
#'
#' @note
#' `rtDuration` considers the total runtime (including MS1 and MS2 scans).
#' 
#' @param spectra `Spectra` object
#'
#' @param msLevel `integer`
#' @param ... not used here
#'
#' @param ... not used here
#'  
#' @return `numeric(4)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}, Johannes Rainer
#' 
#' @export
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847", "unknown"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine", "unknown"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876),
#'     c(83.0603, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994),
#'     c(3.146, 61.611))
#' spd$rtime <- c(9.44, 9.44, 15.84, 15.81)
#' sps <- Spectra(spd)
#' rtOverMsQuarters(spectra = sps, msLevel = 2L)
rtOverMsQuarters <- function(spectra, msLevel = 1L, ...) {

    ## we assume that with RT duration the mzQC consortium means the run time 
    ## of the whole run, including MS1 and MS2
    rtd <- rtDuration(spectra)

    ## truncate spectra based on the msLevel
    spectra <- filterMsLevel(object = spectra, msLevel)
    if (length(spectra) < 4)
        stop("Spectra object does contain less than four spectra")

    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    RT <- rtime(spectra)
    rtmin <- min(RT)
    
    ## partition the spectra (rows) into four parts 
    ## (they are not necessarily equal)
    ind <- sort(rep(seq_len(4), length.out = length(spectra)))
    idx <- which(c(diff(ind), 1) == 1)
    
    ## calculate the retention time inidces
    rt_quarters <- (RT[idx] - rtmin) / rtd
    names(rt_quarters) <- c("Quarter1", "Quarter2", "Quarter3", "Quarter4")
    rt_quarters
}

#' @name ticQuantileToQuantileLogRatio
#' 
#' @title MS1 quantile TIC change ratio to Quantile 1 (QC:4000057) or to 
#' previous quantile (QC:4000058)
#' 
#' @description 
#' "The log ratio for the second to n-th quantile of TIC changes over first 
#' quantile of TIC changes." [PSI:QC]
#' id: QC:4000057
#' 
#' "The log ratio for the second to n-th quantile of TIC over the previous 
#' quantile of TIC. For the boundary elements min/max are used." [PSI:QC]
#' id: QC:4000058
#'
#' @note
#'
#' This function interprets the *quantiles* from the [PSI:QC] definition as
#' *quartiles*, i.e. the 0, 25, 50, 75 and 100% quantiles are used.
#' 
#' @details
#' 
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' is_a: QC:4000023 ! MS1 metric
#'
#' *TIC changes* are interpreted as follows:
#' (1) the cumulative sum (`cumsum`) of the  spectras' TIC is calculated 
#' (with spectra ordered by retention time),
#' 
#' (2) quartiles are then calculated on these, 
#' 
#' (3) for *QC:4000057* the log2 ratio between the 25, 50, 75
#' and 100% quantile to the 0% quantile is calculated. For *QC:4000058* 
#' ratios between the 25/0, 50/25, 75/50 and 100/75% quantiles are calculated.
#' 
#' The `log2` values are returned instead of the `log` values.
#' 
#' @param spectra `Spectra` object
#' @param relativeTo `character(1)`, one of `"Q1"` or `"previous"`
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' sps <- Spectra(spd)
#' ticQuantileToQuantileLogRatio(spectra = sps, relativeTo = "Q1",
#'     msLevel = 2L)
#' ticQuantileToQuantileLogRatio(spectra = sps, relativeTo = "previous",
#'     msLevel = 2L)
ticQuantileToQuantileLogRatio <- function(spectra, relativeTo = "Q1", 
    msLevel = 1L, ...) {
  
    if (length(relativeTo) != 1)
        stop("'relativeTo' has to be of length 1")
    else
        relativeTo <- match.arg(relativeTo, choices = c("Q1", "previous"))
    
    spectra <- filterMsLevel(object = spectra, msLevel)    
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra")
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    
    ## create cumulative sum of tic/ionCount
    TIC <- ionCount(spectra)
    ticSum <- cumsum(TIC)    
    quantileTICSum <- quantile(ticSum, na.rm = TRUE)

    ## calculate the changes in TIC per quantile
    changeQ1 <- quantileTICSum[["25%"]] - quantileTICSum[["0%"]]
    changeQ2 <- quantileTICSum[["50%"]] - quantileTICSum[["25%"]]
    changeQ3 <- quantileTICSum[["75%"]] - quantileTICSum[["50%"]]
    changeQ4 <- quantileTICSum[["100%"]] - quantileTICSum[["75%"]]
    if (relativeTo == "Q1") {
        ratioQuantileTIC <- c(changeQ2, changeQ3, changeQ4) / changeQ1
        names(ratioQuantileTIC) <- c("Q2/Q1", "Q3/Q1", "Q4/Q1")
    }

    ## calculate the ratio between Q2/Q3/Q4 to previous quantile TIC changes
    if (relativeTo == "previous") {
        ratioQuantileTIC <- c(changeQ2 / changeQ1, changeQ3 / changeQ2, 
                            changeQ4 / changeQ3)
        names(ratioQuantileTIC) <- c("Q2/Q1", "Q3/Q2", "Q4/Q3")
    }
    ## take the log2 and return
    log2(ratioQuantileTIC)
}

#' @name numberSpectra
#'
#' @title Number of MS1 or MS2 spectra (QC:4000059/QC:4000060)
#'
#' @description
#' "The number of MS1 events in the run." [PSI:QC]
#' id: QC:4000059
#' "The number of MS2 events in the run." [PSI:QC]
#' id: QC:4000060
#'
#' For *QC:4000059*, `msLevel` is set to 1. For *QC:4000060*, `msLevel` is 
#' set to 2.
#'
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000023 ! MS1 metric
#'
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#'
#' @return `numeric(1)`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @export
#'
#' @importFrom ProtGenerics filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' sps <- Spectra(spd)
#' numberSpectra(spectra = sps, msLevel = 1L)
#' numberSpectra(spectra = sps, msLevel = 2L)
numberSpectra <- function(spectra, msLevel = 1L, ...) {  
    spectra <- filterMsLevel(object = spectra, msLevel)
    length(spectra)
}


#' @name medianPrecursorMz
#' 
#' @title Precursor median m/z for IDs (QC:4000065)
#' 
#' @description
#' "Median m/z value for all identified peptides (unique ions) after 
#' FDR." [PSI:QC]
#' id: QC:4000065
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor m/z values are obtained,
#' 
#' (3) the median value is returned (`NAs` are removed).
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000023 ! MS1 metric
#' is_a: QC:4000025 ! ion source metric
#' 
#' @note
#' `medianPrecursorMz` will calculate the *precursor* median m/z of all 
#' Spectra within `spectra`. If the calculation needs be done according to
#' *QC:4000065*, the `Spectra` object should be prepared accordingly, i.e.
#' filtered with e.g. [filterPrecursorMz()] or subsetted to spectra with
#' identification data.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorMz
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorMz <- c(170.16, 170.16, 195.0876)
#' sps <- Spectra(spd)
#' medianPrecursorMz(spectra = sps, msLevel = 2L)
medianPrecursorMz <- function(spectra, msLevel = 1L, ...) {  
    spectra <- filterMsLevel(object = spectra, msLevel)
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra")
    
    ## FDR correction??
    mz <- precursorMz(spectra)
    median(mz, na.rm = TRUE)
}

#' @name rtIqr
#' 
#' @title Interquartile RT period for peptide identifications (QC:4000072)
#' 
#' @description
#' "The interquartile retention time period, in seconds, for all peptide 
#' identifications over the complete run." [PSI:QC]
#' id: QC:4000072
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the retention time values are obtained,
#' 
#' (3) the interquartile range is obtained from the values and returned
#' (`NA` values are removed).
#'
#' @details
#' Longer times indicate better chromatographic separation.
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' Retention time values that are `NA` are removed.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000072*, the 
#' `Spectra` object should be prepared accordingly, i.e. subsetted to spectra
#' with identification data.
#' 
#' The stored retention time information in `spectra` might have a different
#' unit than seconds. `rtIqr` will return the IQR based on the values stored
#' in `spectra` and will not convert these values to seconds. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics filterMsLevel rtime
#' @importFrom stats IQR
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' rtIqr(spectra = sps, msLevel = 2L)
rtIqr <- function(spectra, msLevel = 1L, ...) {
    spectra <- filterMsLevel(object = spectra, msLevel)    
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra") 
    
    ## get the retention time
    rt <- rtime(spectra)
    
    ## remove the retention time values that are NA and return the interquartile 
    ## range 
    IQR(rt, na.rm = TRUE)
}

#' @name rtIqrRate
#' 
#' @title Peptide identification rate of the interquartile RT period (QC:4000073)
#' 
#' @description
#' "The identification rate of peptides for the interquartile retention time 
#' period, in peptides per second." [PSI:QC]
#' id: QC:4000073
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the retention time values are obtained,
#' 
#' (3) the 25% and 75% quantiles are obtained from the retention time values
#' (`NA` values are removed),
#' 
#' (4) the number of eluted features between this 25% and 75% quantile is 
#' calculated,
#' 
#' (5) the number of features is divided by the interquartile range of the 
#' retention time and returned.
#' 
#' 
#' @details
#' Higher rates indicate efficient sampling and identification.
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000073*, the 
#' `Spectra` object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data.
#' 
#' The stored retention time information in `spectra` might have a different
#' unit than seconds. `rtIqr` will return the IQR based on the values stored
#' in `spectra` and will not convert these values to seconds. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(2)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics rtime
#' @importFrom stats quantile
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' rtIqrRate(spectra = sps, msLevel = 2L)
rtIqrRate <- function(spectra, msLevel = 1L, ...) {
    spectra <- filterMsLevel(object = spectra, msLevel)    
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra") 
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    RT <- rtime(spectra)
    
    quantileRT <- quantile(RT, na.rm = TRUE)
    
    ## get the RT values of the 25% and 75% quantile
    quantile25RT <- quantileRT[["25%"]]
    quantile75RT <- quantileRT[["75%"]]
    
    ## get the number of eluted features between the 25% and 75% quantile
    nFeatures <- RT >= quantile25RT & RT <= quantile75RT
    nFeatures <- sum(nFeatures)
    
    ## divide the number of eluted features between the 25% and 75% quantile
    ## by the IQR to get the elution rate per second 
    nFeatures / rtIqr(spectra, msLevel = msLevel)
}

#' @name areaUnderTic
#' 
#' @title Area under TIC (QC:4000077)
#' 
#' @description 
#' "The area under the total ion chromatogram." [PSI:QC]
#' id: QC:4000077
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the sum of the ion counts are obtained and returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' The sum of the TIC is returned as an equivalent to the area.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics tic ionCount
#' 
#' @examples 
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' sps <- Spectra(spd)
#' areaUnderTic(spectra = sps, msLevel = 2L)
areaUnderTic <- function(spectra, msLevel = 1L, ...) {  
    spectra <- filterMsLevel(object = spectra, msLevel)    
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra") 

    TIC <- ionCount(spectra)
    
    ## sum up the TIC (equivalent to the area) and return
    sum(TIC, na.rm = TRUE)
}

#' @name areaUnderTicRtQuantiles
#' 
#' @title Area under TIC RT quantiles (QC:4000078)
#' 
#' @description 
#' "The area under the total ion chromatogram of the retention time quantiles. 
#' Number of quantiles are given by the n-tuple." [PSI:QC]
#' id: QC:4000078
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the `Spectra` object is ordered according to the retention time,
#' 
#' (3) the 0%, 25%, 50%, 75%, and 100% quantiles of the retention time values
#' are obtained,
#' 
#' (4) the ion count of the intervals between the 0%/25%, 25%/50%, 50%/75%, and
#' 75%/100% are obtained 
#' 
#' (5) the ion counts of the intervals are summed (TIC) and the values returned
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' The sum of the TIC is returned as an equivalent to the area.
#' 
#' @note
#' This function interprets the *quantiles* from the [PSI:QC] definition as
#' *quartiles*, i.e. the 0, 25, 50, 75 and 100% quantiles are used.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(4)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics tic ionCount rtime
#' @importFrom stats quantile
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' areaUnderTicRtQuantiles(spectra = sps, msLevel = 2L)
areaUnderTicRtQuantiles <- function(spectra, msLevel = 1L, ...) {
    spectra <- filterMsLevel(object = spectra, msLevel)
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra")

    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    rt <- rtime(spectra)
    
    quantileRT <- quantile(rt, na.rm = TRUE)
    
    tic <- ionCount(spectra)
    
    ## get the TICs for the 1st, 2nd, 3rd, and 4th quartile
    ticQ1 <- tic[rt > quantileRT[["0%"]] & rt <= quantileRT[["25%"]]]
    ticQ2 <- tic[rt > quantileRT[["25%"]] & rt <= quantileRT[["50%"]]]
    ticQ3 <- tic[rt > quantileRT[["50%"]] & rt <= quantileRT[["75%"]]]
    ticQ4 <- tic[rt > quantileRT[["75%"]] & rt <= quantileRT[["100%"]]]
    
    ## sum the TICs (area) for the 1st, 2nd, 3rd, and 4th quartile
    areaTicQ1 <- sum(ticQ1, na.rm = TRUE)
    areaTicQ2 <- sum(ticQ2, na.rm = TRUE)
    areaTicQ3 <- sum(ticQ3, na.rm = TRUE)
    areaTicQ4 <- sum(ticQ4, na.rm = TRUE)
    
    ## return the summed TICs as a named vector
    areaTic <- c(areaTicQ1, areaTicQ2, areaTicQ3, areaTicQ4)
    names(areaTic) <- c("25%", "50%", "75%", "100%")
    
    areaTic
}

#' @name extentIdentifiedPrecursorIntensity
#' 
#' @title Extent of identified precursor intensity (QC:4000125)
#' 
#' @description 
#' "Ratio of 95th over 5th percentile of precursor intensity for identified 
#' peptides" [PSI:QC]
#' id: QC:4000125
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensities of the precursor ions are obtained,
#' 
#' (3) the 5% and 95% quantile of these intensities are obtained 
#' (`NA` values are removed),
#' 
#' (4) the ratio between the 95% and the 5% intensity quantile is calculated
#' and returned.
#' 
#' @details
#' Can be used to approximate the dynamic range of signal
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' Precursor intensity values that are `NA` are removed.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000125*, the 
#' `Spectra` object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics filterMsLevel precursorIntensity
#' @importFrom stats quantile
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorIntensity <- c(100, 100, 100)
#' sps <- Spectra(spd)
#' extentIdentifiedPrecursorIntensity(spectra = sps, msLevel = 2L)
extentIdentifiedPrecursorIntensity <- function(spectra, msLevel = 1L, ...) {  
    spectra <- filterMsLevel(object = spectra, msLevel)
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra") 
  
    ## retrieve the precursorIntensity and calculate the 5% and 95% quantile
    precInt <- precursorIntensity(spectra)
    quantilePrecInt <- quantile(precInt, probs = c(0.05, 0.95), na.rm = TRUE)
    
    ## calculate the ratio between the 95% and 5% quantile and return the value
    quantilePrecInt[["95%"]] / quantilePrecInt[["5%"]]
}

#' @name medianTicRtIqr
#' 
#' @title Median of TIC values in the RT range in which the middle half of 
#' peptides are identified (QC:4000130)
#' 
#' @description
#' "Median of TIC values in the RT range in which half of peptides are 
#' identified (RT values of Q1 to Q3 of identifications)" [PSI:QC]
#' id: QC:4000130
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the `Spectra` object is ordered according to the retention time,
#' 
#' (3) the features between the 1st and 3rd quartile are obtained 
#' (half of the features that are present in `spectra`),
#' 
#' (4) the ion count of the features within the 1st and 3rd quartile is 
#' obtained,
#' 
#' (5) the median value of the ion count is calculated (`NA` values are 
#' removed) and the median value is returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' The function `medianTicRtIqr` uses the function [ionCount()] as an 
#' equivalent to the TIC.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000130*, the 
#' `Spectra` object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @export
#' 
#' @importFrom ProtGenerics filterMsLevel tic ionCount
#' @importFrom stats median
#'
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' medianTicRtIqr(spectra = sps, msLevel = 2L)
medianTicRtIqr <- function(spectra, msLevel = 1L, ...) {
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0)
        stop("'spectra' does not contain any spectra") 
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    
    ## get the Q1 to Q3 of identifications 
    ## (half of peptides that are identitied)
    ind <- rep(seq_len(4), length.out = length(spectra))
    ind <- sort(ind)
    Q1ToQ3 <- spectra[ind %in% c(2, 3), ]
    
    ## take the ionCount of the Q1 to Q3 of identifications
    ticQ1ToQ3 <- ionCount(Q1ToQ3)
    
    ## take the median value of the TIC within this interval and return it
    median(ticQ1ToQ3, na.rm = TRUE)
}

#' @name medianTicOfRtRange
#' 
#' @title Median of TIC values in the shortest RT range in which half of the 
#' peptides are identified (QC:4000132)
#' 
#' @description 
#' "Median of TIC values in the shortest RT range in which half of the 
#' peptides are identified"  [PSI:QC] 
#' id: QC:4000132
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the `Spectra` object is ordered according to the retention time,
#' 
#' (3) the number of features in `spectra` is obtained and the number for 
#' half of the features is calculated,
#' 
#' (4) iterate through the features (always by taking the neighbouring
#' half of features) and calculate the retention time range of the
#' set of features,
#' 
#' (5) retrieve the set of features with the minimum retention time 
#' range,
#' 
#' (6) calculate from the set of (5) the median TIC (`NA` values are removed)
#' and return it.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' The function `medianTicOfRtRange` uses the function `ionCount` as an 
#' equivalent to the TIC.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000132*, the 
#' `Spectra` object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return
#' `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics filterMsLevel tic ionCount rtime
#' @importFrom stats median
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' medianTicOfRtRange(spectra = sps, msLevel = 2L)
medianTicOfRtRange <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    } 
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)
    rt <- rtime(spectra)
    
    tic <- ionCount(spectra)
    
    ## retrieve number of features in object and calculate the number for half
    ## of the features
    n <- length(spectra)
    n_half <- ceiling(n / 2)
    
    ## iterate through the bunches of features (always take n_half features),
    ## start with 1 + n_half, 2 + n_half, 3 + n_half, ..., i_end + n_half
    ## calculate the RT range
    rangeRT <- lapply(seq_len(n_half + 1), function(i) {
        ind <- seq(i, i - 1 + n_half)
        rt_i <- rt[ind]
        max(rt_i) - min(rt_i)
    })
    rangeRT <- unlist(rangeRT)
    
    ## retrieve the index of the bunch with the minimum range and calculate
    ## the median TIC
    indMin <- which.min(rangeRT)
    ind <- seq(indMin, indMin - 1 + n_half)
    ticMin <- tic[ind]
    median(ticMin, na.rm = TRUE)
}

#' @name mzAcquisitionRange
#' 
#' @title m/z acquisition range (QC:4000138)
#' 
#' @description 
#' "Upper and lower limit of m/z values at which spectra are recorded." [PSI:QC]
#' id: QC:4000138
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the m/z values of the peaks within `spectra` are obtained,
#' 
#' (3) the minimum and maximum m/z values are obtained and returned. 
#'
#' @details
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' is_a: QC:4000004 ! n-tuple
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return
#' `numeric(2)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics filterMsLevel mz
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' sps <- Spectra(spd)
#' mzAcquisitionRange(spectra = sps, msLevel = 2L)
mzAcquisitionRange <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    } 

    mzList <- mz(spectra)
    mz <- unlist(mzList)
    mzRange <- range(mz)
    names(mzRange) <- c("min", "max")
    
    return(mzRange)
}

#' @name rtAcquisitionRange
#' 
#' @title Retention time acquisition range (QC:4000139)
#' 
#' @description 
#' "Upper and lower limit of time at which spectra are recorded." [PSI:QC]
#' id: QC:4000139
#' 
#' #' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the retention time values of the features within `spectra` are obtained,
#' 
#' (3) the minimum and maximum retention time values are obtained and returned. 
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return
#' `numeric(2)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics rtime filterMsLevel
#' 
#' @examples 
#' library(S4Vectors)
#' library(Spectra)
#' 
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' rtAcquisitionRange(spectra = sps, msLevel = 2L)
rtAcquisitionRange <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
    
    rt <- rtime(spectra)
    rtRange <- range(rt)
    names(rtRange) <- c("min", "max")
    
    return(rtRange)
}

#' @name precursorIntensityRange
#' 
#' @title Precursor intensity range (QC:4000144)
#' 
#' @description 
#' "Minimum and maximum precursor intensity recorded." [PSI:QC]
#' id: QC:4000144
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensity of the precursor ions within `spectra` are obtained,
#' 
#' (3) the minimum and maximum precursor intensity values are obtained and 
#' returned. 
#' 
#' @details
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' is_a: QC:4000004 ! n-tuple
#' 
#' The intensity range of the precursors informs about the dynamic range of 
#' the acquisition. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return
#' `numeric(2)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorIntensity filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorIntensity <- c(100.0, 100.0, 100.0)
#' sps <- Spectra(spd)
#' precursorIntensityRange(spectra = sps, msLevel = 2L)
precursorIntensityRange <- function(spectra, msLevel = 1, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
  
    int <- precursorIntensity(spectra)
    rangeInt <- range(int)
    names(rangeInt) <- c("min", "max")
    
    return(rangeInt)
}

#' @name precursorIntensityQuartiles
#' 
#' @title Precursor intensity distribution Q1, Q2, Q3 (QC:4000167), 
#' Identified precursor intensity distribution Q1, Q2, Q3 (QC:4000228), or
#' Unidentified precursor intensity distribution Q1, Q2, Q3 (QC:4000233)
#' 
#' @description
#' "From the distribution of precursor intensities, the quartiles Q1, Q2, Q3" 
#' [PSI:QC]
#' id: QC:4000167
#' 
#' "From the distribution of identified precursor intensities, the quartiles 
#' Q1, Q2, Q3" [PSI:QC]
#' id: QC:40000228
#' 
#' "From the distribution of unidentified precursor intensities, the 
#' quartiles Q1, Q2, Q3" [PSI:QC]
#' id: QC:4000233
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensity of the precursor ions within `spectra` are obtained,
#' 
#' (3) the 25%, 50%, and 75% quantile of the  precursor intensity values are 
#' obtained (`NA` values are removed) and returned.
#' 
#' @details
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' is_a: QC:4000004 ! n-tuple
#' 
#' The intensity distribution of the precursors informs about the dynamic 
#' range of the acquisition.
#' 
#' The intensity distribution of the identified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' The intensity distribution of the unidentified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' @note 
#' The `Spectra` object might contain features that were (not) identified. If
#' the calculation needs to be done according to *QC:4000228*/*QC:4000233*, the 
#' `Spectra` object should be prepared accordingly. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(3)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorIntensity filterMsLevel
#' @importFrom stats quantile
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorIntensity <- c(100.0, 100.0, 100.0)
#' sps <- Spectra(spd)
#' precursorIntensityQuartiles(spectra = sps, msLevel = 2L)
precursorIntensityQuartiles <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
  
    int <- precursorIntensity(spectra)
    quantile(int, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
}


#' @name precursorIntensityMean
#' 
#' @title Precursor intensity distribution mean (QC:4000168),
#' Identified precursor intensity distribution mean (QC:4000229), or
#' Unidentified precursor intensity distribution mean (QC:4000234)
#' 
#' @description
#' "From the distribution of precursor intensities, the mean." [PSI:QC]
#' id: QC:4000168
#' 
#' "From the distribution of identified precursor intensities, the mean"
#'  [PSI:QC]
#' id: QC:4000229
#' 
#' "From the distribution of unidentified precursor intensities, the mean" 
#' [PSI:QC]
#' id: QC:4000234
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensity of the precursor ions within `spectra` are obtained,
#' 
#' (3) the mean of the precursor intensity values is obtained 
#' (`NA` values are removed) and returned. 
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' 
#' The intensity distribution of the precursors informs about the dynamic 
#' range of the acquisition.
#' 
#' The intensity distribution of the identified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' The intensity distribution of the unidentified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' @note 
#' The `Spectra` object might contain features that were (not) identified. If
#' the calculation needs to be done according to *QC:4000229*/*QC:4000234*, the 
#' `Spectra` object should be prepared accordingly. 
#' 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorIntensity filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorIntensity <- c(100.0, 100.0, 100.0)     
#' sps <- Spectra(spd)
#' precursorIntensityMean(spectra = sps, msLevel = 2L)
precursorIntensityMean <- function(spectra, msLevel = 1L, ...) {
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
  
    int <- precursorIntensity(spectra)
    mean(int, na.rm = TRUE)
}

#' @name precursorIntensitySd
#' 
#' @title Precursor intensity distribution sigma (QC:4000169),
#' Identified precursor intensity distribution sigma (QC:4000230), or 
#' Unidentified precursor intensity distribution sigma (QC:4000235)
#' 
#' @description 
#' "From the distribution of precursor intensities, the sigma value." [PSI:QC]
#' id: QC:4000169
#' 
#' "From the distribution of identified precursor intensities, the sigma 
#' value" [PSI:QC]
#' id: QC:4000230
#' 
#' "From the distribution of unidentified precursor intensities, the 
#' sigma value" [PSI:QC]
#' id: QC:4000235
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensity of the precursor ions within `spectra` are obtained,
#' 
#' (3) the standard deviation of precursor intensity values is obtained 
#' (`NA` values are removed) and returned. 
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' 
#' The intensity distribution of the precursors informs about the dynamic 
#' range of the acquisition.
#' 
#' The intensity distribution of the identified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' The intensity distribution of the unidentified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability.
#' 
#' @note 
#' The `Spectra` object might contain features that were (not) identified. If
#' the calculation needs to be done according to *QC:4000230*/*QC:4000235*, the 
#' `Spectra` object should be prepared accordingly. 
#'
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorIntensity filterMsLevel
#' @importFrom stats sd
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorIntensity <- c(100.0, 100.0, 100.0)
#' sps <- Spectra(spd)
#' precursorIntensitySd(spectra = sps, msLevel = 2L)
precursorIntensitySd <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
  
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
    
    int <- precursorIntensity(spectra)
    sd(int, na.rm = TRUE)
}

#' @name msSignal10xChange
#' 
#' @title MS1 signal jump/fall (10x) count (QC:4000172/QC:4000173)
#' 
#' @description 
#' "The count of MS1 signal jump (spectra sum) by a factor of ten or more (10x)
#' between two subsequent scans" [PSI:QC]
#' id: QC:4000172
#' 
#' "The count of MS1 signal decline (spectra sum) by a factor of ten or more 
#' (10x) between two subsequent scans" [PSI:QC]
#' id: QC:4000173
#' 
#' #' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the intensity of the precursor ions within `spectra` are obtained,
#' 
#' (3) the intensity values of the features are obtained via the ion count,
#' 
#' (4) the signal jumps/declines of the intensity values with the two 
#' subsequent intensity values is calculated,
#' 
#' (5) in the case of *QC:4000172*, the signal jumps by a factor of ten or more
#' are counted and returned;
#' in the case of *QC:4000173*, the signal declines by a factor of ten or more
#' are counted and returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' 
#' An unusual high count of signal jumps or falls can indicate ESI stability 
#' issues.
#' 
#' The function `msSignal10xChange` uses the function `ionCount` as an 
#' equivalent to the TIC.
#' 
#' @param spectra `Spectra` object
#' @param change `character(1)`, one of `"jump"` or `"fall"`
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics tic ionCount filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$rtime <- c(9.44, 9.44, 15.84)
#' sps <- Spectra(spd)
#' msSignal10xChange(spectra = sps, change = "jump", msLevel = 2L)
#' msSignal10xChange(spectra = sps, change = "fall", msLevel = 2L)
msSignal10xChange <- function(spectra, change = "jump", msLevel = 1L, ...) {
  
    if (length(change) != 1) {
        stop("'change' has to be of length 1")
    } else {
        change <- match.arg(change, choices = c("jump", "fall"))
    }
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra")
    }
    
    ## order spectra according to increasing retention time
    spectra <- .rtOrderSpectra(spectra)

    tic <- ionCount(spectra)
    
    precedingTic <- tic[seq_len(length(tic) - 1)]
    followingTic <- tic[seq_len(length(tic))[-1]]
    
    ## calculate the ratio between following and preceding TICs and calculate
    ## the number of 10X jumps or falls depending on the change argument
    ratioTic <- followingTic / precedingTic
    
    if (change == "jump")
      numberRatioChange <- sum(ratioTic >= 10)
    if (change == "fall")
      numberRatioChange <- sum(ratioTic <= 0.1)
    
    return(numberRatioChange)
}

#' @name ratioCharge1over2
#' 
#' @title Charged peptides ratio 1+ over 2+ (QC:4000174) or 
#' Charged spectra ratio 1+ over 2+ (QC:4000179)
#' 
#' @description 
#' "Ratio of 1+ peptide count over 2+ peptide count in identified spectra" [PSI:QC]
#' id: QC:4000174
#' 
#' "Ratio of 1+ spectra count over 2+ spectra count in all MS2" [PSI:QC]
#' id: QC:4000179
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor charge is obtained,
#' 
#' (3) the number of precursors with charge 1+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' `NA` is returned if there are no features with precursor charge of 1+ or 
#' 2+.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000174*, the 
#' `Spectra` object should be prepared accordingly. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorCharge <- c(1L, 1L, 1L)
#' sps <- Spectra(spd)
#' ratioCharge1over2(spectra = sps, msLevel = 2L)
ratioCharge1over2 <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(1, 2) %in% names(chargeTable)))
        chargeRatio <- chargeTable[["1"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name ratioCharge3over2
#' 
#' @title Charged peptides ratio 3+ over 2+ (QC:4000175) or 
#' charged spectra ratio 3+ over 2+ (QC:4000180)
#' 
#' @description 
#' "Ratio of 3+ peptide count over 2+ peptide count in identified spectra" [PSI:QC]
#' id: QC:4000175
#' 
#' "Ratio of 3+ peptide count over 2+ peptide count in all MS2" [PSI:QC]
#' id: QC:4000180
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor charge is obtained,
#' 
#' (3) the number of precursors with charge 3+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' `NA` is returned if there are no features with precursor charge of 2+ or 
#' 3+.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000175*, the 
#' `Spectra` object should be prepared accordingly. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorCharge <- c(1L, 1L, 1L)
#' sps <- Spectra(spd)
#' ratioCharge3over2(spectra = sps, msLevel = 2L)
ratioCharge3over2 <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(2, 3) %in% names(chargeTable)))
        chargeRatio <- chargeTable[["3"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name ratioCharge4over2
#' 
#' @title Charged peptides ratio 4+ over 2+ (QC:4000176) or charged spectra 
#' ratio 4+ over 2+ (QC:4000181)
#' 
#' @description 
#' "Ratio of 4+ peptide count  over 2+ peptide count in identified 
#' spectra" [PSI:QC]
#' id: QC:4000176
#' 
#' "Ratio of 4+ peptide count over 2+ peptide count in all MS2" [PSI:QC]
#' id: QC:4000181
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor charge is obtained,
#' 
#' (3) the number of precursors with charge 4+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned.
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000176*, the 
#' `Spectra` object should be prepared accordingly.
#' 
#' `NA` is returned if there are no features with precursor charge of 2+ or 
#' 3+.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorCharge <- c(1L, 1L, 1L)
#' sps <- Spectra(spd)
#' ratioCharge4over2(spectra = sps, msLevel = 2L)
ratioCharge4over2 <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(2, 4) %in% names(chargeTable)))
        chargeRatio <- chargeTable[["4"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}


#' @name meanCharge
#' 
#' @title Mean charge in identified spectra (QC:4000177) or mean precursor 
#' charge in all MS2 (QC:4000182)
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor charge is obtained,
#' 
#' (3) the mean of the precursor charge values is calculated and returned.
#' 
#' @description 
#' "Mean charge in identified spectra" [PSI:QC]
#' id: QC:4000177
#' 
#' "Mean precursor charge in all MS2" [PSI:QC]
#' id: QC:4000182
#'
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000177*, the 
#' `Spectra` object should be prepared accordingly. 
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' spd$precursorCharge <- c(1L, 1L, 1L)
#' sps <- Spectra(spd)
#' meanCharge(spectra = sps, msLevel = 2L)
meanCharge <- function(spectra, msLevel = 1L, ...) {
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' object does not contain any spectra") 
    }
    
    charge <- precursorCharge(spectra)
    mean(charge, na.rm = TRUE)
}

#' @name medianCharge
#' 
#' @title Median charge in identified spectra (QC:4000178) or median
#' precursor charge in all MS2 (QC:4000183)
#' 
#' @description 
#' "Median charge in identified spectra" [PSI:QC]
#' id: QC:4000178
#' 
#' "Median precursor charge in all MS2" [PSI:QC]
#' id: QC:4000183
#' 
#' The metric is calculated as follows:
#' (1) the `Spectra` object is filtered according to the MS level,
#' 
#' (2) the precursor charge is obtained,
#' 
#' (3) the median of the precursor charge values is calculated and returned.
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @note 
#' The `Spectra` object might contain features that were not identified. If
#' the calculation needs to be done according to *QC:4000178*, the 
#' `Spectra` object should be prepared accordingly.
#' 
#' @param spectra `Spectra` object
#' @param msLevel `integer`
#' @param ... not used here
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge filterMsLevel
#' 
#' @examples
#' library(S4Vectors)
#' library(Spectra)
#'
#' spd <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     polarity = c(1L, 1L, 1L),
#'     id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
#'     name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))
#' ## Assign m/z and intensity values
#' spd$mz <- list(
#'     c(109.2, 124.2, 124.5, 170.16, 170.52),
#'     c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
#'     c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
#'         111.0551, 123.0429, 138.0662, 195.0876))
#' spd$intensity <- list(
#'     c(3.407, 47.494, 3.094, 100.0, 13.240),
#'     c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
#'     c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994))
#' sps <- Spectra(spd)
#' spd$precursorCharge <- c(1L, 1L, 1L)
#' medianCharge(spectra = sps, msLevel = 2L)
medianCharge <- function(spectra, msLevel = 1L, ...) {
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("'spectra' does not contain any spectra") 
    }
    
    charge <- precursorCharge(spectra)
    median(charge, na.rm = TRUE)
}

