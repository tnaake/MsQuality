#' @name rtDuration
#' 
#' @title RT duration (QC:4000053)
#' 
#' @description
#' "The retention time duration of the MS run in seconds, similar to the 
#' highest scan time minus the lowest scan time." [PSI:QC]
#' id: QC:4000053
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' 
#' @param spectra `Spectra` object
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
rtDuration <- function(spectra) {
  
    RT <- rtime(object = spectra)
    max(RT) - min(RT)
}

#' @name rtOverTICquantile
#' 
#' @title RT over TIC quantile (QC:4000054)
#' 
#' @description 
#' "The interval when the respective quantile of the TIC accumulates divided by 
#' retention time duration. The number of quantiles observed is given by the 
#' size of the tuple." [PSI:QC]
#' id: QC:4000054
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
#' @param msLevel `numeric(1)`
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
#' rtOverTICquantile(spectra = sps, msLevel = 2L)
rtOverTICquantile <- function(spectra, probs = seq(0, 1, 0.25),
                              msLevel = 1L, ...) {
    
    ## truncate spectra based on the MS level
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    ## order spectra according to increasing retention time
    spectra <- .rt_order_spectra(spectra)
    RT <- rtime(spectra)
    rtmin <- min(RT)
    rtd <- rtDuration(spectra)

    ## My undestanding:
    ## what is the relative duration of the LC run after which the cumulative
    ## TIC exceeds (for the first time) the respective quantile of the
    ## cumulative TIC. The "accumulates" puzzled me a little but I guess you
    ## were right with the assumption they mean the "sum of the TICs of all
    ## previous spectra".
    TIC <- cumsum(ionCount(spectra))
    tq <- quantile(TIC, probs = probs, ...)
    idxs <- unlist(lapply(tq, function(z) which.max(TIC >= z)))
    res <- (RT[idxs] - rtmin) / rtd
    names(res) <- names(idxs)
    res
}

.rt_order_spectra <- function(x) {
    RT <- rtime(x)
    if (!any(is.na(RT)) && is.unsorted(RT))
        x <- x[order(RT)]
    x
}

#' @name rtOverMSQuarters
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
#' @details
#' 
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' is_a: QC:4000023 ! MS1 metric
#'
#' @note
#'
#' RT-Duration considers the total runtime (including MS1 and MS2 scans).
#' 
#' @param spectra `Spectra` object
#'
#' @param msLevel `numeric(1)`
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
#' rtOverMSQuarters(spectra = sps, msLevel = 2L)
rtOverMSQuarters <- function(spectra, msLevel = 1L) {

    ## I would assume that with RT duration they mean the run time of the whole
    ## run, including MS1 and MS2
    rtd <- rtDuration(spectra)
    
    ## truncate spectra based on the MSLevel
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
    
    ## order spectra according to increasing retention time
    spectra <- .rt_order_spectra(spectra)
    RT <- rtime(spectra)
    rtmin <- min(RT)
    
    ## partition the spectra (rows) into four parts 
    ## (they are not necessarily equal)
    ind <- sort(rep(seq_len(4), length.out = length(spectra)))
    idx <- which(c(diff(ind), 1) == 1)
    (RT[idx] - rtmin) / rtd
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
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' is_a: QC:4000023 ! MS1 metric
#'
#' The `log2` values are returned instead of the `log` values.
#'
#' *TIC changes* are interpreted as the cumulative sum (`cumsum`) of the
#' spectras' TIC (with spectra ordered by retention time). Quartiles are then
#' calculated on these. For *QC:4000057* the log2 ratio between the 25, 50, 75
#' and 100% quartile to the 0% quartile is calculated. For *QC:4000058* 
#' ratios between the 25/0, 50/25, 75/50 and 100/75% quartiles are calculated.
#' 
#' @param spectra `Spectra` object
#' 
#' @param relativeTo `character`
#'
#' @param msLevel `integer(1)`
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
#'     MSLevel = 2L)
#' ticQuantileToQuantileLogRatio(spectra = sps, relativeTo = "previous",
#'     MSLevel = 2L)
ticQuantileToQuantileLogRatio <- function(spectra, 
                                          relativeTo = c("Q1", "previous"),
                                          msLevel = 1L) {
  
    relativeTo <- match.arg(relativeTo)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## order spectra according to increasing retention time
    spectra <- .rt_order_spectra(spectra)
    
    ## create cumulative sum of tic/ionCount
    TIC <- cumsum(ionCount(spectra))
    
    quantileTICSum <- quantile(TIC, na.rm = TRUE)

    if (relativeTo == "Q1") {
        ## Not sure what they mean by "TIC changes", thus calculating here
        ## only ratio between TIC.
        ratioQuantileTIC <- quantileTICSum[c("25%", "50%", "75%", "100%")] /
            quantileTICSum["0%"]
        names(ratioQuantileTIC) <- c("Q1/Q0", "Q2/Q0", "Q3/Q0", "Q4/Q0")
    } else {
        ratioQuantileTIC <- quantileTICSum[c("25%", "50%", "75%", "100%")] /
            quantileTICSum[c("0%", "25%", "50%", "75%")]
        names(ratioQuantileTIC) <- c("Q1/Q0", "Q2/Q1", "Q3/Q2", "Q4/Q3")
    }
    log2(ratioQuantileTIC)
}

#' @name numberSpectra
#' 
#' @title Number of MS1/MS2 spectra (QC:4000059/QC:4000060)
#' 
#' @description
#' "The number of MS1 events in the run." [PSI:QC]
#' id: QC:4000059
#' "The number of MS2 events in the run." [PSI:QC]
#' id: QC:4000060
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000023 ! MS1 metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric(1)`
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
#' numberSpectra(spectra = sps, MSLevel = 1L)
#' numberSpectra(spectra = sps, MSLevel = 2L)
numberSpectra <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    len <- length(spectra)
    
    return(len)
}


#' @name medianPrecursorMZ
#' 
#' @title Precursor median m/z for IDs (QC:4000065)
#' 
#' @description
#' "Median m/z value for all identified peptides (unique ions) after 
#' FDR." [PSI:QC]
#' id: QC:4000065
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000023 ! MS1 metric
#' is_a: QC:4000025 ! ion source metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' medianPrecursorMZ(spectra = sps, MSLevel = 2L)
medianPrecursorMZ <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ################ FDR correction??????????? 
    mz <- ProtGenerics::precursorMz(spectra)
    medianMZ <- median(mz, na.rm = TRUE)
    
    return(medianMZ)
}

#' @name rtIQR
#' 
#' @title Interquartile RT period for peptide identifications (QC:4000072)
#' 
#' @description
#' "The interquartile retention time period, in seconds, for all peptide 
#' identifications over the complete run." [PSI:QC]
#' id: QC:4000072
#'
#' @details
#' Longer times indicate better chromatographic separation.
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' rtIQR(spectra = sps, MSLevel = 2L)
rtIQR <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## get the retention time
    RT <- ProtGenerics::rtime(spectra)
    
    ## IQR???, what is the unit for rtime, always seconds??
    iqr <- stats::IQR(RT, na.rm = TRUE)
    
    return(iqr)
}

#' @name rtIQRrate
#' 
#' @title Peptide identification rate of the interquartile RT period (QC:4000073)
#' 
#' @description
#' "The identification rate of peptides for the interquartile retention time 
#' period, in peptides per second." [PSI:QC]
#' id: QC:4000073
#' 
#' @details
#' Higher rates indicate efficient sampling and identification.
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' rtIQRrate(spectra = sps, MSLevel = 2L)
rtIQRrate <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## order spectra according to increasing retention time
    RT <- ProtGenerics::rtime(spectra)
    spectra <- spectra[order(RT)]
    RT <- RT[order(RT)]
    
    quantileRT <- stats::quantile(RT, na.rm = TRUE)
    
    ## get the RT values of the 25% and 75% quantile
    quantile25RT <- quantileRT[["25%"]]
    quantile75RT <- quantileRT[["75%"]]
    
    ## get the number of eluted features between the 25% and 75% quantile
    nFeatures <- RT >= quantile25RT & RT <= quantile75RT
    nFeatures <- sum(nFeatures)
    
    ## divide the number of eluted features between the 25% and 75% quantile
    ## by the IQR to get the elution rate per second 
    rate <- nFeatures / rtIQR(spectra, MSLevel = MSLevel)
    
    return(rate)
}

#' @name areaUnderTIC
#' 
#' @title Area under TIC (QC:4000077)
#' 
#' @description 
#' "The area under the total ion chromatogram." [PSI:QC]
#' id: QC:4000077
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' The sum of the TIC is returned as an equivalent to the area.
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' areaUnderTIC(spectra = sps, MSLevel = 2L)
areaUnderTIC <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }

    TIC <- ProtGenerics::ionCount(spectra)
    
    ## sum up the TIC (equivalent to the area) and return
    area <- sum(TIC, na.rm = TRUE)
    
    return(area)
}

#' @name areaUnderTICRTquantiles
#' 
#' @title Area under TIC RT quantiles (QC:4000078)
#' 
#' @description 
#' "The area under the total ion chromatogram of the retention time quantiles. 
#' Number of quantiles are given by the n-tuple." [PSI:QC]
#' id: QC:4000078
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' The sum of the TIC is returned as an equivalent to the area.
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' areaUnderTICRTquantiles(spectra = sps, MSLevel = 2L)
areaUnderTICRTquantiles <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }

    ## order spectra according to increasing retention time
    RT <- ProtGenerics::rtime(spectra)
    spectra <- spectra[order(RT)]
    RT <- RT[order(RT)]
    
    quantileRT <- stats::quantile(RT, na.rm = TRUE)
    
    TIC <- ProtGenerics::ionCount(spectra)
    
    ## get the TICs for the 1st, 2nd, 3rd, and 4th quartile
    ticQ1 <- TIC[RT > quantileRT[["0%"]] & RT <= quantileRT[["25%"]]]
    ticQ2 <- TIC[RT > quantileRT[["25%"]] & RT <= quantileRT[["50%"]]]
    ticQ3 <- TIC[RT > quantileRT[["50%"]] & RT <= quantileRT[["75%"]]]
    ticQ4 <- TIC[RT > quantileRT[["75%"]] & RT <= quantileRT[["100%"]]]
    
    ## sum the TICs (area) for the 1st, 2nd, 3rd, and 4th quartile
    areaTICQ1 <- sum(ticQ1, na.rm = TRUE)
    areaTICQ2 <- sum(ticQ2, na.rm = TRUE)
    areaTICQ3 <- sum(ticQ3, na.rm = TRUE)
    areaTICQ4 <- sum(ticQ4, na.rm = TRUE)
    
    ## return the summed TICs as a named vector
    areaTIC <- c(areaTICQ1, areaTICQ2, areaTICQ3, areaTICQ4)
    names(areaTIC) <- c("25%", "50%", "75%", "100%")
    
    return(areaTIC)
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
#' @details
#' Can be used to approximate the dynamic range of signal
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' extentIdentifiedPrecursorIntensity(spectra = sps, MSLevel = 2L)
extentIdentifiedPrecursorIntensity <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
  
    ## retrieve the precursorIntensity and calculate the 5% and 95% quantile
    precInt <- ProtGenerics::precursorIntensity(spectra)
    quantilePrecInt <- stats::quantile(precInt, probs = c(0.05, 0.95), 
                                                              na.rm = TRUE)
    
    ## calculate the ratio between the 95% and 5% quantile and return the value
    ratio <- quantilePrecInt[["95%"]] / quantilePrecInt[["5%"]]
    
    return(ratio)
}

#' @name medianTICRTIQR
#' 
#' @title Median of TIC values in the RT range in which the middle half of 
#' peptides are identified (QC:4000130)
#' 
#' @description
#' "Median of TIC values in the RT range in which half of peptides are 
#' identified (RT values of Q1 to Q3 of identifications)" [PSI:QC]
#' id: QC:4000130
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' The function `medianTICrtIQR` uses the function `ionCount` as an 
#' equivalent to the TIC.
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' medianTICRTIQR(spectra = sps, MSLevel = 2L)
medianTICRTIQR <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    } 
    ## order spectra according to increasing retention time
    RT <- ProtGenerics::rtime(spectra)
    spectra <- spectra[order(RT)]
    RT <- RT[order(RT)]
    
    ## get the Q1 to Q3 of identifications 
    ## (half of peptides that are identitied)
    ind <- rep(seq_len(4), length.out = length(spectra))
    ind <- sort(ind)
    Q1ToQ3 <- spectra[ind %in% c(2, 3), ]
    
    # ## take the TIC of the Q1 to Q3 of identifications
    
    ## how to define TIC? take the ionCount?? #######################
    ## take the ionCount of the Q1 to Q3 of identifications
    ticQ1ToQ3 <- ProtGenerics::ionCount(Q1ToQ3)
    
    ## take the median value of the TIC within this interval and return
    medianTICQ1ToQ3 <- stats::median(ticQ1ToQ3)
    
    return(medianTICQ1ToQ3)
}

#' @name medianTICofRTRange
#' 
#' @title Median of TIC values in the shortest RT range in which half of the 
#' peptides are identified (QC:4000132)
#' 
#' @description 
#' "Median of TIC values in the shortest RT range in which half of the 
#' peptides are identified"  [PSI:QC] 
#' id: QC:4000132
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric

#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' medianTICofRTRange(spectra = sps, MSLevel = 2L)
medianTICofRTRange <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    } 
    
    ## order spectra according to increasing retention time
    RT <- ProtGenerics::rtime(spectra)
    spectra <- spectra[order(RT)]
    RT <- RT[order(RT)]
    
    TIC <- ProtGenerics::ionCount(spectra)
    
    ## retrieve number of features in object and calculate the number for half
    ## of the features
    n <- length(spectra)
    n_half <- ceiling(n / 2)
    
    ## iterate through the bunches of features (always take n_half features),
    ## start with 1 + n_half, 2 + n_half, 3 + n_half, ..., i_end + n_half
    ## calculate the RT range
    rangeRT <- lapply(seq_len(n_half + 1), function(i) {
        ind <- seq(i, i - 1 + n_half)
        rt <- RT[ind]
        max(rt) - min(rt)
    })
    rangeRT <- unlist(rangeRT)
    
    ## retrieve the index of the bunch with the minimum range and calculate
    ## the median TIC
    indMin <- which.min(rangeRT)
    ind <- seq(indMin, indMin - 1 + n_half)
    ticMin <- TIC[ind]
    medianTICMin <- stats::median(ticMin, na.rm = TRUE)
    
    return(medianTICMin)
}

#' @name mzAcquisitionRange
#' 
#' @title m/z acquisition range (QC:4000138)
#' 
#' @description 
#' "Upper and lower limit of m/z values at which spectra are recorded." [PSI:QC]
#' id: QC:4000138
#' 
#' @details
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' is_a: QC:4000004 ! n-tuple
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' mzAcquisitionRange(spectra = sps, MSLevel = 2L)
mzAcquisitionRange <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    } 
    
    mzList <- ProtGenerics::mz(spectra)
    MZ <- unlist(mzList)
    mzRange <- range(MZ)
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
#' @details
#' is_a: QC:4000004 ! n-tuple
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' rtAcquisitionRange(spectra = sps, MSLevel = 2L)
rtAcquisitionRange <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
    
    RT <- ProtGenerics::rtime(spectra)
    rtRange <- range(RT)
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
#' @details
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' is_a: QC:4000004 ! n-tuple
#' 
#' The intensity range of the precursors informs about the dynamic range of 
#' the acquisition. 
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' precursorIntensityRange(spectra = sps, MSLevel = 2L)
precursorIntensityRange <- function(spectra, MSLevel = 1) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
  
    int <- ProtGenerics::precursorIntensity(spectra)
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
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' precursorIntensityQuartiles(spectra = sps, MSLevel = 2L)
precursorIntensityQuartiles <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
  
    int <- ProtGenerics::precursorIntensity(spectra)
    quartiles <- stats::quantile(int, probs = c(0.25, 0.50, 0.75), 
                                                                na.rm = TRUE)
    
    return(quartiles)
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
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' precursorIntensityMean(spectra = sps, MSLevel = 2L)
precursorIntensityMean <- function(spectra, MSLevel = 1L) {
    
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
  
    int <- ProtGenerics::precursorIntensity(spectra)
    intMean <- mean(int, na.rm = TRUE)
    
    return(intMean)
}

#' @name precursorIntensitySD
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
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' precursorIntensitySD(spectra = sps, MSLevel = 2L)
precursorIntensitySD <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
  
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
    
    int <- ProtGenerics::precursorIntensity(spectra)
    mzSD <- stats::sd(int, na.rm = TRUE)
    
    return(mzSD)
}

#' @name msSignal10XChange
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
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' An unusual high count of signal jumps or falls can indicate ESI stability 
#' issues.
#' 
#' @param spectra `Spectra` object
#' @param change `character`, one of `"jump"` or `"fall"`
#' @param MSLevel `numeric`
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
#' sps <- Spectra(spd)
#' msSignal10XChange(spectra = sps, change = "jump", MSLevel = 2L)
#' msSignal10XChange(spectra = sps, change = "fall", MSLevel = 2L)
msSignal10XChange <- function(spectra, change = c("jump", "fall"), 
                                                          MSLevel = 1L) {

    change <- match.arg(change)
    
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra")
    }
    
    ## order spectra according to increasing retention time
    RT <- ProtGenerics::rtime(spectra)
    spectra <- spectra[order(RT)]
    RT <- RT[order(RT)]

    TIC <- ProtGenerics::ionCount(spectra)
    
    precedingTIC <- TIC[seq_len(length(TIC) - 1)]
    followingTIC <- TIC[seq_len(length(TIC))[-1]]
    
    ## calculate the ratio between following and preceding TICs and calculate
    ## the number of 10X jumps or falls depending on the change argument
    ratioTIC <- followingTIC / precedingTIC
    
    if (change == "jump")
      numberRatioChange <- sum(ratioTIC >= 10)
    if (change == "fall")
      numberRatioChange <- sum(ratioTIC <= 0.1)
    
    return(numberRatioChange)
}

#' @name RatioCharge1over2
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
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' RatioCharge1over2(spectra = sps, MSLevel = 2L)
RatioCharge1over2 <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- ProtGenerics::precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(1, 2) %in% names(chargeTable)))
        chargeRatio <- chargeTable[["1"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name RatioCharge3over2
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
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' RatioCharge3over2(spectra = sps, MSLevel = 2L)
RatioCharge3over2 <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- ProtGenerics::precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(2, 3) %in% names(chargeTable)))
        chargeRatio <- chargeTable[["3"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name RatioCharge4over2
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
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' RatioCharge4over2(spectra = sps, MSLevel = 2L)
RatioCharge4over2 <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    ## is there a way to get charge of actual entries, not only of precursor?
    charge <- ProtGenerics::precursorCharge(spectra)
    
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
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' meanCharge(spectra = sps, MSLevel = 2L)
meanCharge <- function(spectra, MSLevel = 1L) {
    
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    charge <- ProtGenerics::precursorCharge(spectra)
    chargeMean <- mean(charge, na.rm = TRUE)
    return(chargeMean)
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
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' @param MSLevel `numeric`
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
#' medianCharge(spectra = sps, MSLevel = 2L)
medianCharge <- function(spectra, MSLevel = 1L) {
  
    spectra <- ProtGenerics::filterMsLevel(object = spectra, MSLevel)
    
    if (length(spectra) == 0) {
        stop("Spectra object does not contain any spectra") 
    }
    
    charge <- ProtGenerics::precursorCharge(spectra)
    chargeMedian <- median(charge, na.rm = TRUE)
    return(chargeMedian)
}

