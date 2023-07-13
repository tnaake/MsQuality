#' @name chromatographyDuration
#' 
#' @title chromatography duration (MS:4000053)
#' 
#' @description
#' MS:4000053
#' "The retention time duration of the chromatography in seconds." [PSI:MS] \cr 
#' 
#' The metric is calculated as follows: \cr 
#' (1) the retention time associated to the individual Spectra is obtained, \cr  
#' (2) the maximum and the minimum of the retention time is obtained, \cr  
#' (3) the difference between the maximum and the minimum is calculated and 
#' returned. \cr 
#' 
#' @details
#' MS:4000053
#' synonym: "RT-Duration" RELATED [PMID:24494671] \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000016 ! retention time metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_value_concept NCIT:C25330 ! Duration \cr
#' relationship: has_units UO:0000010 ! second \cr
#' 
#' Retention time values that are \code{NA} are removed.
#' 
#' @param spectra \code{Spectra} object
#' 
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
#' chromatographyDuration(spectra = sps)
chromatographyDuration <- function(spectra, ...) {  
    RT <- rtime(object = spectra)
    res <- max(RT, na.rm = TRUE) - min(RT, na.rm = TRUE)
    
    ## add attributes and return
    attributes(res) <- list(chromatographyDuration = "MS:4000053")
    res
}

#' @name ticQuartersRtFraction
#' 
#' @title TIC quarters RT fraction (MS:4000054)
#' 
#' @description
#' MS:4000054
#' "The interval when the respective quarter of the TIC accumulates divided by 
#' retention time duration." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is ordered according to the retention time, \cr 
#' (2) the cumulative sum of the ion count is calculated (TIC), \cr 
#' (3) the quantiles are calculated according to the \code{probs} argument, e.g.
#' when \code{probs} is set to \code{c(0, 0.25, 0.5, 0.75, 1)} the 0\%, 25\%, 50\%, 75\% 
#' and 100\% quantile is calculated, \cr 
#' (4) the retention time/relative retention time (retention time divided by 
#' the total run time taking into account the minimum retention time) is 
#' calculated, \cr 
#' (5) the (relative) duration of the LC run after which the cumulative
#' TIC exceeds (for the first time) the respective quantile of the
#' cumulative TIC is calculated and returned. \cr
#' 
#' @details
#' MS:4000054
#' synonym: "RT-TIC-Q1" RELATED [PMID:24494671] \cr
#' synonym: "RT-TIC-Q2" RELATED [PMID:24494671] \cr
#' synonym: "RT-TIC-Q3" RELATED [PMID:24494671] \cr
#' synonym: "RT-TIC-Q4" RELATED [PMID:24494671] \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000016 ! retention time metric \cr
#' relationship: has_metric_category MS:4000017 ! chromatogram metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units UO:0000191 ! fraction \cr
#' 
#' @param spectra \code{Spectra} object
#' @param probs \code{numeric} defining the quantiles. See \code{probs = seq(0, 1, 0.25)}.
#' @param msLevel \code{integer}
#' @param relative \code{logical}, if set to \code{TRUE} the relative retention time 
#' will be returned instead of the abolute retention time
#' @param ... not used here
#' 
#' @return \code{numeric} of length equal to length \code{probs} with the relative
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
#' ticQuartersRtFraction(spectra = sps, msLevel = 2L)
ticQuartersRtFraction <- function(spectra, probs = seq(0, 1, 0.25),
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
    ## probs * max(TIC) calculates the portion of the total TIC at the probs,
    ## e.g. if sum of TIC is 10000 and probs = c(0, 1, 0.25), ids will be 
    ## the indices where the TIC is first higher than 0, 2500, 5000, 7500, 
    ## and 10000
    TIC <- cumsum(ionCount(spectra))
    idxs <- lapply(probs * max(TIC), function(z) which(TIC >= z)[1]) |>
        unlist()
    
    ##idxs <- unlist(lapply(ticQuantile, function(z) which.max(TIC >= z)))
    
    if (relative) {
        rtMin <- min(RT)
        chromatographyDuration <- chromatographyDuration(spectra)
        res <- (RT[idxs] - rtMin) / chromatographyDuration  
    } else {
        res <- RT[idxs]
    }
    
    ## add attributes and return
    attributes(res) <- list(names = paste0(probs * 100, "%"), 
        ticQuartersRtFraction = "MS:4000054")
    res
}

#' @title Order Spectra according to increasing retention time
#' 
#' @description 
#' The function \code{.rtOrderSpectra} orders the features in a \code{Spectra} object
#' according to the (increasing) retention time values. 
#' 
#' @details
#' Internal function in quality metric functions.
#' 
#' @param spectra \code{Spectra} object
#' 
#' @return \code{Spectra} object with the features ordered according to the 
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
#' MsQuality:::.rtOrderSpectra(sps)
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
#' @title MS1 quarter RT fraction (MS:4000055) or 
#' MS2 quarter RT fraction (MS:4000056)
#' 
#' @description
#' MS:4000055 \cr
#' "The interval used for acquisition of the first, second, third, and fourth 
#' quarter of all MS1 events divided by retention time duration." [PSI:MS] \cr
#' 
#' MS:4000056 \cr
#' "The interval used for acquisition of the first, second, third, and fourth 
#' quarter of all MS2 events divided by retention time duration." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the retention time duration of the whole \code{Spectra} object is determined
#' (taking into account all the MS levels), \cr 
#' (2) the \code{Spectra} object is filtered according to the MS level and 
#' subsequently ordered according to the retention time \cr 
#' (3) the MS events are split into four (approximately) equal parts, \cr 
#' (4) the relative retention time is calculated (using the retention time 
#' duration from (1) and taking into account the minimum retention time), \cr 
#' (5) the relative retention time values associated to the MS event parts
#' are returned.
#' 
#' @details
#' MS:4000055 \cr
#' synonym: "RT-MS-Q1" RELATED [PMID:24494671] \cr
#' synonym: "RT-MS-Q2" RELATED [PMID:24494671] \cr
#' synonym: "RT-MS-Q3" RELATED [PMID:24494671] \cr
#' synonym: "RT-MS-Q4" RELATED [PMID:24494671] \cr
#' is_a: MS:4000004 ! n-tuple \cr 
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000016 ! retention time metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units UO:0000191 ! fraction \cr
#'
#' MS:4000056 \cr
#' synonym: "RT-MSMS-Q1" RELATED [PMID:24494671] \cr
#' synonym: "RT-MSMS-Q2" RELATED [PMID:24494671] \cr
#' synonym: "RT-MSMS-Q3" RELATED [PMID:24494671] \cr
#' synonym: "RT-MSMS-Q4" RELATED [PMID:24494671] \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000016 ! retention time metric \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units UO:0000191 ! fraction \cr
#' 
#' The function returns \code{c(NaN, NaN, NaN, NaN)} if the filtered 
#' \code{spectra} object has less than 4 scan events.
#'
#' An attribute will only be returned if \code{msLevel} is 1 or 2.
#' 
#' @note
#' \code{chromatographyDuration} considers the total runtime (including MS1 
#' and MS2 scans).
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' @param ... not used here
#'  
#' @return \code{numeric(4)}
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

    if (length(msLevel) != 1)
        stop("'msLevel' has to be of length 1")
    
    ## we assume that with RT duration the mzQC consortium means the run time 
    ## of the whole run, including MS1 and MS2
    rtd <- chromatographyDuration(spectra)

    ## truncate spectra based on the msLevel
    spectra <- filterMsLevel(object = spectra, msLevel)
    if (length(spectra) < 4) {
        res <- c(NaN, NaN, NaN, NaN)
    } else {
        ## order spectra according to increasing retention time
        spectra <- .rtOrderSpectra(spectra)
        RT <- rtime(spectra)
        rtmin <- min(RT)
        
        ## partition the spectra (rows) into four parts 
        ## (they are not necessarily equal)
        ind <- sort(rep(seq_len(4), length.out = length(spectra)))
        idx <- which(c(diff(ind), 1) == 1)
        
        ## calculate the retention time inidces
        res <- (RT[idx] - rtmin) / rtd
    }
    
    names(res) <- c("Quarter1", "Quarter2", "Quarter3", "Quarter4")
    
    ## add attributes and return
    if (msLevel == 1L) 
        ms_term <- "MS:4000055"
    if (msLevel == 2L) 
        ms_term <- "MS:4000056"
    if (msLevel %in% c(1L, 2L))
        attributes(res) <- c(attributes(res), list(rtOverMsQuarters = ms_term))
    
    res
}

#' @name ticQuartileToQuartileLogRatio
#' 
#' @title MS1 TIC-change quartile ratios (MS:4000057) or  
#' MS1 TIC quartile ratios (MS:4000058)
#' 
#' @description 
#' MS:4000057 \cr
#' "The log ratios of successive TIC-change quartiles. The TIC changes are 
#' the list of MS1 total ion current (TIC) value changes from one to the next 
#' scan, produced when each MS1 TIC is subtracted from the preceding MS1 TIC. 
#' The metric's value triplet represents the log ratio of the TIC-change 
#' Q2 to Q1, Q3 to Q2, TIC-change-max to Q3" [PSI:MS] \cr 
#' For calculation of MS:400057 set \code{mode = "TIC_change"}. \cr
#' 
#' MS:4000058 \cr
#' The log ratios of successive TIC quartiles. The metric's value triplet 
#' represents the log ratios of TIC-Q2 to TIC-Q1, TIC-Q3 to TIC-Q2, 
#' TIC-max to TIC-Q3." [PSI:MS] \cr
#' For calculation of MS:400058 set \code{mode = "TIC"}. \cr
#'
#' The metric is calculated as follows: \cr
#' (1) the TIC (\code{ionCount}) of the  \code{spectra} is calculated per scan
#' event (with spectra ordered by retention time), \cr
#' (2) for *MS:4000057*, the differences between TIC values are calculated 
#' between subsequent scan events,
#' for *MS:4000058*, the TIC values between subsequent scan events are taken
#' as they are, \cr
#' (3) for *MS:4000057* and *MS:4000058* the ratios between the 25\%, 50\%, 
#' 75\%, and 100\% quantile to the 25% quantile of the values of (2) are 
#' calculated.
#' Alternatively, if \code{relativeTo = "Q1"}, the ratios are calculated 
#' between the 50\%/25\%, 75\%/25\%, and 100\%/25\% quantiles, \cr
#' (4) The \code{log} values of the ratios are returned. \cr
#' 
#' @note
#' This function interprets the *quantiles* from the [PSI:QC] definition as
#' *quartiles*, i.e. the 0, 25, 50, 75 and 100\% quantiles are used.
#' 
#' @details
#' MS:4000057 \cr
#' synonym: "MS1-TIC-Change-Q2" RELATED [PMID:24494671] \cr
#' synonym: "MS1-TIC-Change-Q3" RELATED [PMID:24494671] \cr
#' synonym: "MS1-TIC-Change-Q4" RELATED [PMID:24494671] \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000017 ! chromatogram metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_value_concept STATO:0000105 ! log signal intensity ratio \cr
#' 
#' MS:4000058 \cr
#' synonym: "MS1-TIC-Q2" RELATED [PMID:24494671] \cr
#' synonym: "MS1-TIC-Q3" RELATED [PMID:24494671] \cr
#' synonym: "MS1-TIC-Q4" RELATED [PMID:24494671] \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000017 ! chromatogram metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_value_concept STATO:0000105 ! log signal intensity ratio \cr
#' 
#' An attribute will only be returned if \code{relativeTo} is 
#' \code{"previous"} and \code{msLevel} is 1. 
#'  
#' @param spectra \code{Spectra} object
#' @param relativeTo \code{character(1)}, one of \code{"Q1"} or 
#' \code{"previous"}
#' @param mode \code{character(1)}, one of \code{"TIC_change"} or \code{"TIC"}
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
#' 
#' ## MS:4000057
#' ticQuartileToQuartileLogRatio(spectra = sps, relativeTo = "previous",
#'     msLevel = 2L, mode = "TIC_change")
#' ticQuartileToQuartileLogRatio(spectra = sps, relativeTo = "Q1",
#'     msLevel = 2L, mode = "TIC_change")
#' 
#' ## MS:4000058
#' ticQuartileToQuartileLogRatio(spectra = sps, relativeTo = "previous",
#'     msLevel = 2L, mode = "TIC")
#' ticQuartileToQuartileLogRatio(spectra = sps, relativeTo = "Q1",
#'     msLevel = 2L, mode = "TIC")
ticQuartileToQuartileLogRatio <- function(spectra, 
    relativeTo = c("previous", "Q1"), mode = "TIC_change", msLevel = 1L, ...) {

    if (length(relativeTo) != 1)
        stop("'relativeTo' has to be of length 1")
    else
        relativeTo <- match.arg(relativeTo)
    
    if (length(mode) != 1)
        stop("'relativeTo' has to be of length 1")
    else
        mode <- match.arg(mode, choices = c("TIC_change", "TIC"))
    
    spectra <- filterMsLevel(object = spectra, msLevel)    
    if (length(spectra) == 0) {
        
        ratioQuantileTIC <- c(NaN, NaN, NaN)
        
    } else {
        ## order spectra according to increasing retention time
        spectra <- .rtOrderSpectra(spectra)
        
        ## calculate TIC/ionCount
        TIC <- ionCount(spectra)
        
        ## calculate TIC changes if mode == "TIC_change",
        ## otherwise if mode == "TIC" take values as they are
        if (mode == "TIC_change")
            TIC <- diff(TIC) 
        if (mode == "TIC")
            TIC <- TIC
        
        ## obtain the quantiles (will be quartiles, 0%, 25%, 50%, 75%, 100%)
        quantileTIC <- quantile(TIC, na.rm = TRUE)
        
        ## calculate the changes in TIC per quantile
        changeQ1 <- quantileTIC[["25%"]]
        changeQ2 <- quantileTIC[["50%"]] 
        changeQ3 <- quantileTIC[["75%"]]
        changeQ4 <- quantileTIC[["100%"]]
        
        ## calculate the ratio between Q2/Q3/Q4 to Q1 quantile TIC (changes)
        if (relativeTo == "Q1") {
            ratioQuantileTIC <- c(changeQ2, changeQ3, changeQ4) / changeQ1
        }
        
        ## calculate the ratio between Q2/Q3/Q4 to previous quantile TIC (changes)
        if (relativeTo == "previous") {
            ratioQuantileTIC <- c(changeQ2 / changeQ1, changeQ3 / changeQ2, 
                                  changeQ4 / changeQ3)
        }
    }
    
    ## add names 
    if (relativeTo == "Q1") {
        names(ratioQuantileTIC) <- c("Q2/Q1", "Q3/Q1", "Q4/Q1")
    }
    if (relativeTo == "previous") {
        names(ratioQuantileTIC) <- c("Q2/Q1", "Q3/Q2", "Q4/Q3")
    }

    ## take the log and return
    res <- log(ratioQuantileTIC)
    
    ## add attributes and return
    if (mode == "TIC_change" & msLevel == 1L) 
        ms_term <- "MS:4000057"
    if (mode == "TIC" & msLevel == 1L) 
        ms_term <- "MS:4000058"
    
    ## add the attribute only in case of relativeTo = "previous"
    if (relativeTo == "previous" & msLevel == 1L)
        attributes(res) <- c(attributes(res), 
            list(ticQuartileToQuartileLogRatio = ms_term))
    
    res
}

#' @name numberSpectra
#'
#' @title number of MS1 spectra (MS:4000059) or number of MS2 spectra 
#' (MS:4000060)
#'
#' @description
#' MS:4000059 \cr
#' "The number of MS1 events in the run." [PSI:MS] \cr
#' 
#' MS:4000060 \cr
#' "The number of MS2 events in the run." [PSI:MS] \cr
#'
#' For *MS:4000059*, \code{msLevel} is set to 1. For *MS:4000060*, 
#' \code{msLevel} is set to 2.
#'
#' @details
#' MS:4000059 \cr
#' synonym: "MS1-Count" EXACT [PMID:24494671] \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_value_type xsd:int ! The allowed value-type for this CV term \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' 
#' MS:4000060 \cr
#' synonym: "MS2-Count" EXACT [PMID:24494671] \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_value_type xsd:int ! The allowed value-type for this CV term \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' 
#' An attribute will only be returned if \code{msLevel} is either 1 or 2.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#'
#' @return \code{numeric(1)}
#'
#' @author Thomas Naake
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
    
    if (length(msLevel) != 1)
        stop("'msLevel' has to be of length 1")
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    res <- length(spectra)
    
    ## add attributes and return
    if (msLevel == 1L) 
        ms_term <- "MS:4000059"
    if (msLevel == 2L) 
        ms_term <- "MS:4000060"
    if (msLevel %in% c(1L, 2L))
        attributes(res) <- list(numberSpectra = ms_term)
    
    res
}

#' @name medianPrecursorMz
#' 
#' @title precursor median m/z for IDs (MS:4000152)
#' 
#' @description
#' MS:4000152 \cr
#' "Median m/z value for all identified peptides (unique ions) after 
#' FDR." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the precursor m/z values are obtained, \cr 
#' (3) the median value is returned (\code{NA}s are removed). \cr
#' 
#' @details
#' MS:4000152 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000021 ! MS1 metric \cr
#' is_a: MS:4000020 ! ion source metric \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is
#' \code{"identified"} and \code{msLevel} is 1.
#' 
#' @note
#' \code{medianPrecursorMz} will calculate the *precursor* median m/z of all 
#' Spectra within \code{spectra}. If the calculation needs be done according to
#' *MS:4000152*, the \code{Spectra} object should be prepared accordingly, i.e.
#' filtered with e.g. [filterPrecursorMz()] or subsetted to spectra with
#' identification data.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
#' 
#' @author Thomas Naake
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
medianPrecursorMz <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        mz <- precursorMz(spectra)
        res <- median(mz, na.rm = TRUE)    
    }
    
    ## add attributes and return
    if (identificationLevel == "identified" & msLevel == 1L)
        attributes(res) <- list(medianPrecursorMz = "MS:4000152")
    
    res
    
}

#' @name rtIqr
#' 
#' @title interquartile retention time period for peptide identifications 
#' (MS:4000153)
#' 
#' @description
#' MS:4000153 \cr
#' "The interquartile retention time period, in seconds, for all peptide 
#' identifications over the complete run. Longer times indicate better 
#' chromatographic separation." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr
#' (2) the retention time values are obtained, \cr
#' (3) the interquartile range is obtained from the values and returned
#' (\code{NA} values are removed). \cr
#'
#' @details
#' MS:4000153 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000017 ! chromatogram metric \cr
#' 
#' Retention time values that are \code{NA} are removed.
#' 
#' An attribute will only be returned if \code{identificationLevel} is
#' \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000153*, the 
#' \code{Spectra} object should be prepared accordingly, i.e. subsetted to 
#' \code{spectra} with identification data.
#' 
#' The stored retention time information in \code{spectra} might have a 
#' different unit than seconds. \code{rtIqr} will return the IQR based on the 
#' values stored in \code{spectra} and will not convert these values to seconds. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
rtIqr <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel) 
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        ## get the retention time
        rt <- rtime(spectra)
        
        ## remove the retention time values that are NA and return the interquartile 
        ## range 
        res <- IQR(rt, na.rm = TRUE)    
    }
    
    ## add attributes and return
    if (identificationLevel == "identified")
        attributes(res) <- list(rtIqr = "MS:4000153")
    
    res
}

#' @name rtIqrRate
#' 
#' @title peptide identification rate of the interquartile RT period (MS:4000154)
#' 
#' @description
#' MS:4000154 \cr
#' "The identification rate of peptides for the interquartile retention time 
#' period, in peptides per second. Higher rates indicate efficient sampling 
#' and identification." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the retention time values are obtained, \cr 
#' (3) the 25\% and 75\% quantiles are obtained from the retention time values
#' (\code{NA} values are removed), \cr 
#' (4) the number of eluted features between this 25\% and 75\% quantile is 
#' calculated, \cr 
#' (5) the number of features is divided by the interquartile range of the 
#' retention time and returned. \cr 
#' 
#' @details
#' MS:4000154 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000017 ! chromatogram metric \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is 
#' \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000154*, the 
#' \code{Spectra} object should be prepared accordingly, i.e. being subsetted to
#' \code{spectra} with identification data.
#' 
#' The stored retention time information in \code{spectra} might have a different
#' unit than seconds. \code{.rtIqr} will return the IQR based on the values stored
#' in \code{spectra} and will not convert these values to seconds. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(2)}
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
rtIqrRate <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
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
        res <- nFeatures / rtIqr(spectra, msLevel = msLevel)    
    }
    
    ## add attributes and return
    if (identificationLevel == "identified")
        attributes(res) <- list(rtIqrRate = "MS:4000154")
    
    res
}

#' @name areaUnderTic
#' 
#' @title area under TIC (MS:4000155)
#' 
#' @description 
#' MS:4000155 \cr
#' "The area under the total ion chromatogram." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the sum of the ion counts are obtained and returned. 
#' 
#' @details
#' MS:4000155 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:40000MS ! ID free \cr
#' is_a: MS:4000017 ! chromatogram metric \cr
#' 
#' The sum of the TIC is returned as an equivalent to the area. \cr
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        TIC <- ionCount(spectra)
        
        ## sum up the TIC (equivalent to the area) and return
        res <- sum(TIC, na.rm = TRUE)
    }
    
    ## add attributes and return
    attributes(res) <- list(areaUnderTic = "MS:4000155")
    
    res
}

#' @name areaUnderTicRtQuantiles
#' 
#' @title area under TIC RT quantiles (MS:4000156)
#' 
#' @description 
#' MS:4000156 \cr
#' "The area under the total ion chromatogram of the retention time quantiles. 
#' Number of quantiles are given by the n-tuple." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the \code{Spectra} object is ordered according to the retention time, \cr 
#' (3) the 0\%, 25\%, 50\%, 75\%, and 100\% quantiles of the retention time 
#' values are obtained, \cr 
#' (4) the ion count of the intervals between the 0\%/25\%, 25\%/50\%, 
#' 50\%/75\%, and 75\%/100\% are obtained, \cr 
#' (5) the ion counts of the intervals are summed (TIC) and the values returned.
#' 
#' @details
#' MS:4000156 \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' is_a: MS:4000009 ! ID free \cr
#' is_a: MS:4000017 ! chromatogram metric \cr
#' 
#' The sum of the TIC is returned as an equivalent to the area. \cr
#' 
#' @note
#' This function interprets the *quantiles* from the [PSI:QC] definition as
#' *quartiles*, i.e. the 0, 25, 50, 75 and 100\% quantiles are used.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return \code{.numeric(4)}
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
    
    if (length(spectra) == 0) {
        res <- c(NaN, NaN, NaN, NaN)
    } else {
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
        res <- c(areaTicQ1, areaTicQ2, areaTicQ3, areaTicQ4)   
    }
    
    ## add attributes and return
    attributes(res) <- list(
        names = names(res) <- c("25%", "50%", "75%", "100%"),
        areaUnderTicRtQuantiles = "MS:4000156")
    
    res
}

#' @name extentIdentifiedPrecursorIntensity
#' 
#' @title extent of identified precursor intensity (MS:4000157)
#' 
#' @description 
#' MS:4000157 \cr
#' "Ratio of 95th over 5th percentile of precursor intensity for identified 
#' peptides. Can be used to approximate the dynamic range of signal." 
#' [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensities of the precursor ions are obtained, \cr 
#' (3) the 5\% and 95\% quantile of these intensities are obtained 
#' (\code{NA} values are removed), \cr 
#' (4) the ratio between the 95\% and the 5\% intensity quantile is calculated
#' and returned. \cr
#' 
#' @details
#' MS:4000157 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' Precursor intensity values that are \code{NA} are removed. \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is
#' \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000157*, the 
#' \code{Spectra} object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
extentIdentifiedPrecursorIntensity <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) { 
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        ## retrieve the precursorIntensity and calculate the 5% and 95% 
        ## quantile
        precInt <- precursorIntensity(spectra)
        quantilePrecInt <- quantile(precInt, probs = c(0.05, 0.95), 
            na.rm = TRUE)
        
        ## calculate the ratio between the 95% and 5% quantile and return the 
        ## value
        res <- quantilePrecInt[["95%"]] / quantilePrecInt[["5%"]]
    }
    
    ## add attributes and return
    if (identificationLevel == "identified")
        attributes(res) <- list(extentIdentifiedPrecursorIntensity = "MS:4000157")
    
    res
}

#' @name medianTicRtIqr
#' 
#' @title median of TIC values in the RT range in which the middle half of 
#' peptides are identified (MS:4000158)
#' 
#' @description
#' MS:4000158 \cr
#' "Median of TIC values in the RT range in which half of peptides are 
#' identified (RT values of Q1 to Q3 of identifications)" [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the \code{Spectra} object is ordered according to the retention time, \cr 
#' (3) the features between the 1st and 3rd quartile are obtained 
#' (half of the features that are present in \code{spectra}), \cr 
#' (4) the ion count of the features within the 1st and 3rd quartile is 
#' obtained, \cr 
#' (5) the median value of the ion count is calculated (\code{NA} values are 
#' removed) and the median value is returned. \cr
#' 
#' @details
#' MS:4000158 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' The function \code{medianTicRtIqr} uses the function [ionCount()] as an 
#' equivalent to the TIC.
#' 
#' An attribute will only be returned if \code{identificationLevel} is 
#' \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000158*, the 
#' \code{Spectra} object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
medianTicRtIqr <- function(spectra, msLevel = 1L,
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
    
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
        res <- median(ticQ1ToQ3, na.rm = TRUE)
    }
    
    ## add attributes and return
    if (identificationLevel == "identified")
        attributes(res) <- list(medianTicRtIqr = "MS:4000158")
    
    res
}

#' @name medianTicOfRtRange
#' 
#' @title median of TIC values in the shortest RT range in which half of the 
#' peptides are identified (MS:4000159)
#' 
#' @description 
#' MS:4000159 \cr
#' "Median of TIC values in the shortest RT range in which half of the 
#' peptides are identified"  [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the \code{Spectra} object is ordered according to the retention time, \cr 
#' (3) the number of features in \code{spectra} is obtained and the number for 
#' half of the features is calculated, \cr 
#' (4) iterate through the features (always by taking the neighbouring
#' half of features) and calculate the retention time range of the
#' set of features, \cr 
#' (5) retrieve the set of features with the minimum retention time 
#' range, \cr 
#' (6) calculate from the set of (5) the median TIC (\code{NA} values are removed)
#' and return it.
#' 
#' @details
#' MS:4000159 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' The function \code{medianTicOfRtRange} uses the function \code{ionCount} as an 
#' equivalent to the TIC.
#' 
#' An attribute will only be returned if \code{identificationLevel} is 
#' \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000159*, the 
#' \code{Spectra} object should be prepared accordingly, i.e. being subsetted to
#' spectra with identification data. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return
#' \code{numeric(1)}
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
medianTicOfRtRange <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
  
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
        
    } else {
    
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
        res <- median(ticMin, na.rm = TRUE)
    }
    
    ## add attributes and return
    if (identificationLevel == "identified")
        attributes(res) <- list(medianTicOfRtRange = "MS:4000159")
    
    res
}

#' @name mzAcquisitionRange
#' 
#' @title m/z acquisition range (MS:4000069)
#' 
#' @description 
#' MS:4000069 \cr
#' "Upper and lower limit of m/z precursor values at which MSn spectra are 
#' recorded." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the m/z values of the peaks within \code{spectra} are obtained, \cr 
#' (3) the minimum and maximum m/z values are obtained and returned. 
#'
#' @details
#' MS:4000069 \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000019 ! MS metric \cr
#' relationship: has_units MS:1000040 ! m/z \cr
#' relationship: has_value_concept STATO:0000035 ! range \cr
#'  
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return \code{numeric(2)}
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
        res <- c(NaN, NaN)
    } else {
        mzList <- mz(spectra)
        mz <- unlist(mzList)
        res <- range(mz)
    }
    
    ## add attributes and return
    attributes(res) <- list(names = c("min", "max"), 
        mzAcquisitionRange = "MS:4000069")
    res
}


#' @name rtAcquisitionRange
#' 
#' @title retention time acquisition range (MS:4000070)
#' 
#' @description 
#' MS:4000070 \cr
#' "Upper and lower limit of retention time at which spectra are recorded." 
#' [PSI:MS] \cr
#' 
#' #' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the retention time values of the features within \code{spectra} are 
#' obtained, \cr 
#' (3) the minimum and maximum retention time values are obtained and 
#' returned. \cr
#' 
#' @details
#' MS:4000070 \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000016 ! retention time metric \cr
#' relationship: has_units UO:0000010 ! second \cr
#' relationship: has_value_concept STATO:0000035 ! range \cr
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return
#' \code{numeric(2)}
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
        res <- c(NaN, NaN)
    } else {
        rt <- rtime(spectra)
        res <- range(rt)
    }
    
    ## add attributes and return
    attributes(res) <- list(names = c("min", "max"), 
        rtAcquisitionRange = "MS:4000070")
    res
}

#' @name precursorIntensityRange
#' 
#' @title precursor intensity range (MS:4000160)
#' 
#' @description 
#' MS:4000160 \cr
#' "Minimum and maximum precursor intensity recorded." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensity of the precursor ions within \code{spectra} are obtained, \cr 
#' (3) the minimum and maximum precursor intensity values are obtained and 
#' returned. 
#' 
#' @details
#' MS:4000160 \cr
#' is_a: MS:4000009 ! ID free \cr
#' is_a: MS:4000001 ! QC metric \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' 
#' The intensity range of the precursors informs about the dynamic range of 
#' the acquisition. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return
#' \code{numeric(2)}
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
        res <- c(NaN, NaN)
    } else {
        int <- precursorIntensity(spectra)
        res <- range(int)
    }

    ## add names
    
    ## add attributes and return
    attributes(res) <- list(names = c("min", "max"), 
        precursorIntensityRange = "MS:4000160")
    
    res
}

#' @name precursorIntensityQuartiles
#' 
#' @title Precursor intensity distribution Q1, Q2, Q3 (MS:4000116), 
#' Identified precursor intensity distribution Q1, Q2, Q3 (MS:4000161), or
#' Unidentified precursor intensity distribution Q1, Q2, Q3 (MS:4000162)
#' 
#' @description
#' MS:4000116 \cr
#' "From the distribution of precursor intensities, the quantiles. I.e. a value 
#' triplet represents the quartiles Q1, Q2, Q3. The intensity distribution of 
#' the precursors informs about the dynamic range of the acquisition." 
#' [PSI:MS] \cr
#' 
#' MS:40000161 \cr
#' "From the distribution of identified precursor intensities, the quartiles 
#' Q1, Q2, Q3. The intensity distribution of the precursors informs about 
#' the dynamic range of the acquisition in relation to identifiability." 
#' [PSI:MS] \cr
#' 
#' id: MS:4000162 \cr
#' "From the distribution of unidentified precursor intensities, the quartiles 
#' Q1, Q2, Q3. The intensity distribution of the precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability." 
#' [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensity of the precursor ions within \code{spectra} are obtained, \cr 
#' (3) the 25\%, 50\%, and 75\% quantile of the  precursor intensity values are 
#' obtained (\code{NA} values are removed) and returned. \cr
#' 
#' @details
#' id: MS:4000116 \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_value_concept STATO:0000291 ! quantile \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr
#' 
#' MS:4000161 \cr
#' is_a: MS:4000004 ! n-tuple \cr
#' is_a: MS:4000008 ! ID based \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_value_concept STATO:0000291 ! quantile \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr
#' 
#' id: MS:4000162 \cr
#' is_a: MS:4000008 ! ID based \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_value_concept STATO:0000291 ! quantile \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr
#'  
#' @note 
#' The \code{Spectra} object might contain features that were (not) identified. If
#' the calculation needs to be done according to *MS:4000161*/*MS:4000162*, the 
#' \code{Spectra} object should be prepared accordingly. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(3)}
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
#' 
#' precursorIntensityQuartiles(spectra = sps, msLevel = 2L)
precursorIntensityQuartiles <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
  
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- c(NaN, NaN, NaN)
    } else {
        int <- precursorIntensity(spectra)
        res <- quantile(int, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <- "MS:4000116"
    if (identificationLevel == "identified")
        ms_term <- "MS:4000161"
    if (identificationLevel == "unidentified")
        ms_term <- "MS:4000162"
    
    attributes(res) <- list(names = c("Q1", "Q2", "Q3"), 
        precursorIntensityQuartiles = ms_term)
    res
}


#' @name precursorIntensityMean
#' 
#' @title precursor intensity distribution mean (MS:4000117),
#' identified precursor intensity distribution mean (MS:4000163), or
#' unidentified precursor intensity distribution mean (MS:4000164)
#' 
#' @description
#' MS:4000117 \cr
#' "From the distribution of precursor intensities, the mean. The intensity 
#' distribution of the precursors informs about the dynamic range of the 
#' acquisition." [PSI:MS] \cr
#' 
#' MS:4000163 \cr
#' "From the distribution of identified precursor intensities, the mean. The 
#' intensity distribution of the identified precursors informs about the dynamic 
#' range of the acquisition in relation to identifiability." [PSI:MS] \cr
#' 
#' MS:4000164 \cr
#' "From the distribution of unidentified precursor intensities, the mean. The 
#' intensity distribution of the unidentified precursors informs about the 
#' dynamic range of the acquisition in relation to identifiability." 
#' [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensity of the precursor ions within \code{spectra} are obtained, \cr 
#' (3) the mean of the precursor intensity values is obtained 
#' (\code{NA} values are removed) and returned. \cr
#' 
#' @details
#' MS:4000117 \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_value_concept STATO:0000401 ! sample mean \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr
#' 
#' MS:4000163 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' relationship: has_value_concept STATO:0000401 ! sample mean \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr
#' 
#' MS:4000164 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' relationship: has_value_concept STATO:0000401 ! sample mean \cr
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV term \cr
#' relationship: has_units MS:1000043 ! intensity unit \cr 
#' 
#' @note 
#' The \code{Spectra} object might contain features that were (not) identified. If
#' the calculation needs to be done according to *MS:4000163*/*MS:4000164*, the 
#' \code{Spectra} object should be prepared accordingly. 
#' 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
precursorIntensityMean <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        int <- precursorIntensity(spectra)
        res <- mean(int, na.rm = TRUE)
    }
  
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <- "MS:4000117"
    if (identificationLevel == "identified")
        ms_term <- "MS:4000163"
    if (identificationLevel == "unidentified")
        ms_term <- "MS:4000164"
    
    attributes(res) <- list(precursorIntensityMean = ms_term)
    
    res
}

#' @name precursorIntensitySd
#' 
#' @title precursor intensity distribution sigma (MS:4000118),
#' identified precursor intensity distribution sigma (MS:4000165), or 
#' unidentified precursor intensity distribution sigma (MS:4000166)
#' 
#' @description 
#' MS:4000118 \cr
#' "From the distribution of precursor intensities, the sigma value. The 
#' intensity distribution of the precursors informs about the dynamic range 
#' of the acquisition." [PSI:MS] \cr
#' 
#' MS:4000165 \cr
#' "From the distribution of identified precursor intensities, the sigma value. 
#' The intensity distribution of the precursors informs about the dynamic range 
#' of the acquisition in relation to identifiability." [PSI:MS] \cr
#' 
#' MS:4000166 \cr
#' "From the distribution of unidentified precursor intensities, the sigma value. 
#' The intensity distribution of the precursors informs about the dynamic range 
#' of the acquisition in relation to identifiability." [PSI:MS] \cr
#' 
#' 
#' The metric is calculated as follows:
#' 
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensity of the precursor ions within \code{spectra} are obtained, \cr 
#' (3) the standard deviation of precursor intensity values is obtained 
#' (\code{NA} values are removed) and returned. 
#' 
#' @details
#' MS:4000118 \cr 
#' is_a: MS:4000003 ! single value \cr 
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr 
#' relationship: has_value_concept STATO:0000237 ! standard deviation \cr 
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr 
#' relationship: has_units MS:1000043 ! intensity unit \cr 
#' 
#' MS:4000165 \cr 
#' is_a: MS:4000003 ! single value \cr 
#' relationship: has_metric_category MS:4000008 ! ID based \cr 
#' relationship: has_value_concept STATO:0000237 ! standard deviation \cr 
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr 
#' relationship: has_units MS:1000043 ! intensity unit \cr 
#' 
#' MS:4000166 \cr 
#' is_a: MS:4000003 ! single value \cr 
#' relationship: has_metric_category MS:4000008 ! ID based \cr 
#' relationship: has_value_concept STATO:0000237 ! standard deviation \cr 
#' relationship: has_value_type xsd:float ! The allowed value-type for this CV 
#' term \cr 
#' relationship: has_units MS:1000043 ! intensity unit \cr 
#'  
#' @note 
#' The \code{Spectra} object might contain features that were (not) identified. If
#' the calculation needs to be done according to *MS:4000165*/*MS:4000166*, the 
#' \code{Spectra} object should be prepared accordingly. 
#'
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
precursorIntensitySd <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
  
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        int <- precursorIntensity(spectra)
        res <- sd(int, na.rm = TRUE)
    }

    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <- "MS:4000118"
    if (identificationLevel == "identified")
        ms_term <- "MS:4000165"
    if (identificationLevel == "unidentified")
        ms_term <- "MS:4000166"
    
    attributes(res) <- list(precursorIntensitySd = ms_term)
    
    res
}

#' @name msSignal10xChange
#' 
#' @title MS1 signal jump/fall (10x) count (MS:4000097/MS:4000098)
#' 
#' @description 
#' MS:4000097 \cr
#' "The number of times where MS1 TIC increased more than 10-fold between 
#' adjacent MS1 scans. An unusual high count of signal jumps or falls can 
#' indicate ESI stability issues." [PSI:MS] \cr
#' 
#' MS:4000098 \cr
#' "The number of times where MS1 TIC decreased more than 10-fold between 
#' adjacent MS1 scans. An unusual high count of signal jumps or falls can 
#' indicate ESI stability issues." [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the intensity of the precursor ions within \code{spectra} are obtained, \cr 
#' (3) the intensity values of the features are obtained via the ion count, \cr 
#' (4) the signal jumps/declines of the intensity values with the two 
#' subsequent intensity values is calculated, \cr 
#' (5) in the case of *MS:4000097*, the signal jumps by a factor of ten or more
#' are counted and returned; \cr
#' in the case of *MS:4000098*, the signal declines by a factor of ten or more
#' are counted and returned. \cr
#' 
#' @details
#' MS:4000097 \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' relationship: has_value_type xsd:integer ! The allowed value-type for this CV 
#' term \cr
#' synonym: "IS-1A"  RELATED [] \cr
#'  
#' The function \code{msSignal10xChange} uses the function \code{ionCount} as an 
#' equivalent to the TIC.
#' 
#' An attribute will only be returned if \code{msLevel} is 1.
#' 
#' @param spectra \code{Spectra} object
#' @param change \code{character(1)}, one of \code{"jump"} or \code{"fall"}
#' @param msLevel \code{integer}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
        res <- NaN
    } else {
        ## order spectra according to increasing retention time
        spectra <- .rtOrderSpectra(spectra)
        
        tic <- ionCount(spectra)
        
        precedingTic <- tic[seq_len(length(tic) - 1)]
        followingTic <- tic[seq_len(length(tic))[-1]]
        
        ## calculate the ratio between following and preceding TICs and calculate
        ## the number of 10X jumps or falls depending on the change argument
        ratioTic <- followingTic / precedingTic
        
        if (change == "jump") 
            res <- sum(ratioTic >= 10)
        if (change == "fall")
            res <- sum(ratioTic <= 0.1)
    }
    
    ## add attributes and return
    if (change == "jump" & msLevel == 1L)
        ms_term <- "MS:4000097"
    if (change == "fall" & msLevel == 1L)
        ms_term <- "MS:4000098"
    if (msLevel == 1L)
        attributes(res) <- list(msSignal10xChange = ms_term)
    
    res
}

#' @name numberEmptyScans
#'
#' @title number of empty MS1 scans (MS:4000099), number of empty MS2 scans 
#' (MS:4000100), or number of empty MS3 scans (MS:4000101)
#'
#' @description
#' MS:4000099 \cr
#' "Number of MS1 scans where the scans' peaks intensity sums to 0 
#' (i.e. no peaks or only 0-intensity peaks)." [PSI:MS] \cr
#' 
#' MS:4000100 \cr
#' "Number of MS2 scans where the scans' peaks intensity sums to 0 
#' (i.e. no peaks or only 0-intensity peaks)." [PSI:MS] \cr
#' 
#' MS:4000101 \cr
#' "Number of MS3 scans where the scans' peaks intensity sums to 0 
#' (i.e. no peaks or only 0-intensity peaks)." [PSI:MS] \cr
#'
#' @details 
#' MS:4000099 \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000021 ! MS1 metric \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' relationship: has_value_type xsd:integer ! The allowed value-type for this 
#' CV term \cr
#' 
#' MS:4000100 \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_metric_category MS:4000022 ! MS2 metric \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' relationship: has_value_type xsd:integer ! The allowed value-type for this 
#' CV term \cr
#' 
#' MS:4000101 \cr
#' is_a: MS:4000003 ! single value \cr
#' relationship: has_metric_category MS:4000009 ! ID free metric \cr
#' relationship: has_metric_category MS:4000012 ! single run based metric \cr
#' relationship: has_units UO:0000189 ! count unit \cr
#' relationship: has_value_type xsd:integer ! The allowed value-type for this CV 
#' term \cr
#' 
#' #' For *MS:4000099*, \code{msLevel} is set to 1. For *MS:4000100*, \code{msLevel} is 
#' set to 2. For *MS:4000101*, \code{msLevel} is set to 3.
#' 
#' An attribute will only be returned if \code{msLevel} is either 1, 2, or 3.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param ... not used here
#'
#' @return \code{numeric(1)}
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
numberEmptyScans <- function(spectra, msLevel = 1L, ...) {
    
    if (length(msLevel) != 1)
        stop("'msLevel' has to be of length 1")
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    ## three cases to take into account: 1) entry is NULL, 2) entry is NA,
    ## or 3) entry is of length 0; in all three cases set to TRUE, otherwise 
    ## to FALSE
    res <- intensity(spectra) |>
        lapply(FUN = function(i) 
            ifelse(is.null(i), TRUE, is.na(i) | length(i) == 0)) |> 
        unlist() |>
        sum()
    
    ## add attributes and return
    if (msLevel == 1L) 
        ms_term <- "MS:4000099"
    if (msLevel == 2L) 
        ms_term <- "MS:4000100"
    if (msLevel == 3L) 
        ms_term <- "MS:4000101"
    if (msLevel %in% c(1L, 2L, 3L))
        attributes(res) <- list(numberEmptyScans = ms_term)
    
    res
}

#' @name ratioCharge1over2
#' 
#' @title charged peptides ratio 1+ over 2+ (MS:4000167) or 
#' charged spectra ratio 1+ over 2+ (MS:4000168)
#' 
#' @description 
#' MS:4000167 \cr
#' "Ratio of 1+ peptide count over 2+ peptide count in identified spectra" 
#' [PSI:MS] \cr
#' 
#' MS:4000168 \cr
#' "Ratio of 1+ spectra count over 2+ spectra count in all MS2" [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the precursor charge is obtained, \cr 
#' (3) the number of precursors with charge 1+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned. \cr
#' 
#' @details
#' MS:4000167 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' MS:4000168 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000009 ! ID free metric \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' \code{NA} is returned if there are no features with precursor charge of 1+ or 
#' 2+. \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is either 
#' \code{"all"} or \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000167*, the 
#' \code{Spectra} object should be prepared accordingly. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
ratioCharge1over2 <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
  
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        ## is there a way to get charge of actual entries, not only of precursor?
        charge <- precursorCharge(spectra)
        
        ## get the number of precursor per charge
        chargeTable <- table(charge)
        
        if (all(c(1, 2) %in% names(chargeTable)))
            res <- chargeTable[["1"]] / chargeTable[["2"]]
        else 
            res <- NaN
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <-  "MS:4000168"
    if (identificationLevel == "identified")
        ms_term <-  "MS:4000167"
    
    if (identificationLevel %in% c("all", "identified"))
        attributes(res) <- list(ratioCharge1over2 = ms_term)
    
    res
    
}

#' @name ratioCharge3over2
#' 
#' @title charged peptides ratio 3+ over 2+ (MS:4000169) or 
#' charged spectra ratio 3+ over 2+ (MS:4000170)
#' 
#' @description 
#' MS:4000169 \cr
#' "Ratio of 3+ peptide count over 2+ peptide count in identified spectra" 
#' [PSI:QC] \cr
#' 
#' MS:4000170 \cr
#' "Ratio of 3+ peptide count over 2+ peptide count in all MS2" [PSI:QC] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the precursor charge is obtained, \cr 
#' (3) the number of precursors with charge 3+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned. \cr
#' 
#' @details
#' MS:4000169 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' MS:4000170 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000009 ! ID free metric \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' \code{NA} is returned if there are no features with precursor charge of 2+ or 
#' 3+.
#' 
#' An attribute will only be returned if \code{identificationLevel} is either 
#' \code{"all"} or \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000169*, the 
#' \code{Spectra} object should be prepared accordingly. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
ratioCharge3over2 <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
  
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        ## is there a way to get charge of actual entries, not only of precursor?
        charge <- precursorCharge(spectra)
        
        ## get the number of precursor per charge
        chargeTable <- table(charge)
        
        if (all(c(2, 3) %in% names(chargeTable)))
            res <- chargeTable[["3"]] / chargeTable[["2"]]
        else 
            res <- NaN
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <-  "MS:4000170"
    if (identificationLevel == "identified")
        ms_term <-  "MS:4000169"
    if (identificationLevel %in% c("all", "identified"))
        attributes(res) <- list(ratioCharge3over2 = ms_term)
    
    res
}

#' @name ratioCharge4over2
#' 
#' @title charged peptides ratio 4+ over 2+ (MS:4000171) or charged spectra 
#' ratio 4+ over 2+ (MS:4000172)
#' 
#' @description 
#' MS:4000171 \cr
#' "Ratio of 4+ peptide count  over 2+ peptide count in identified 
#' spectra" [PSI:MS] \cr
#' 
#' MS:4000172 \cr
#' "Ratio of 4+ peptide count over 2+ peptide count in all MS2" [PSI:MS] \cr
#'  
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the precursor charge is obtained, \cr 
#' (3) the number of precursors with charge 4+ is divided by the number of 
#' precursors with charge 2+ and the ratio is returned.
#' 
#' @details
#' MS:4000171 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' MS:4000172 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000009 ! ID free metric \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is either 
#' \code{"all"} or \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000171*, the 
#' \code{Spectra} object should be prepared accordingly. \cr
#' 
#' \code{NA} is returned if there are no features with precursor charge of 2+ or 
#' 3+.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
ratioCharge4over2 <- function(spectra, msLevel = 1L,
        identificationLevel = c("all", "identified"), ...) {
  
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        ## is there a way to get charge of actual entries, not only of precursor?
        charge <- precursorCharge(spectra)
        
        ## get the number of precursor per charge
        chargeTable <- table(charge)
        
        if (all(c(2, 4) %in% names(chargeTable)))
            res <- chargeTable[["4"]] / chargeTable[["2"]]
        else 
            res <- NaN
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <-  "MS:4000172"
    if (identificationLevel == "identified")
        ms_term <-  "MS:4000171"
    if (identificationLevel %in% c("all", "identified"))
        attributes(res) <- list(ratioCharge4over2 = ms_term)
    
    res
}


#' @name meanCharge
#' 
#' @title Mean precursor charge in identified spectra (MS:4000173) or mean 
#' precursor charge in all MS2 (MS:4000174)
#' 
#' @description 
#' MS:4000173 \cr
#' "Mean precursor charge in identified spectra" [PSI:MS] \cr
#' 
#' MS:4000174 \cr
#' "Mean precursor charge in all MS2" [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr
#' (2) the precursor charge is obtained, \cr
#' (3) the mean of the precursor charge values is calculated and returned. \cr
#' 
#'
#' @details
#' MS:4000173 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' MS:4000174 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID free metric \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' An attribute will only be returned if \code{identificationLevel} is either 
#' \code{"all"} or \code{"identified"}.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000173*, the 
#' \code{Spectra} object should be prepared accordingly. 
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
meanCharge <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {
    
    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        charge <- precursorCharge(spectra)
        res <- mean(charge, na.rm = TRUE)   
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        ms_term <- "MS:4000174"
    if (identificationLevel == "identified")
        ms_term <- "MS:4000173"
    if (identificationLevel %in% c("all", "identified"))
        attributes(res) <- list(meanCharge = ms_term)
    
    res
    
}

#' @name medianCharge
#' 
#' @title median precursor charge in identified spectra (MS:4000175) or median
#' precursor charge in all MS2 (MS:4000176)
#' 
#' @description 
#' MS:4000175 \cr
#' "Median precursor charge in identified spectra" [PSI:MS] \cr
#' 
#' MS:4000176 \cr
#' "Median precursor charge in all MS2" [PSI:MS] \cr
#' 
#' The metric is calculated as follows: \cr
#' (1) the \code{Spectra} object is filtered according to the MS level, \cr 
#' (2) the precursor charge is obtained, \cr 
#' (3) the median of the precursor charge values is calculated and returned.
#' 
#' @details
#' MS:4000175 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000008 ! ID based \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' MS:4000176 \cr
#' is_a: MS:4000003 ! single value \cr
#' is_a: MS:4000009 ! ID free metric \cr
#' is_a: MS:4000001 ! QC metric \cr
#' 
#' An attribute will only be returned if \code{msLevel} is either 1, 2, or 3.
#' 
#' @note 
#' The \code{Spectra} object might contain features that were not identified. If
#' the calculation needs to be done according to *MS:4000175*, the 
#' \code{Spectra} object should be prepared accordingly.
#' 
#' @param spectra \code{Spectra} object
#' @param msLevel \code{integer}
#' @param identificationLevel \code{character(1)}, one of \code{"all"}, 
#' \code{"identified"}, or \code{"unidentified"}
#' @param ... not used here
#' 
#' @return \code{numeric(1)}
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
medianCharge <- function(spectra, msLevel = 1L, 
        identificationLevel = c("all", "identified", "unidentified"), ...) {

    identificationLevel <- match.arg(identificationLevel)
    
    spectra <- filterMsLevel(object = spectra, msLevel)
    
    if (length(spectra) == 0) {
        res <- NaN
    } else {
        charge <- precursorCharge(spectra)
        res <- median(charge, na.rm = TRUE)     
    }
    
    ## add attributes and return
    if (identificationLevel == "all")
        attributes(res) <- "MS:4000176"
    if (identificationLevel == "identified")
        attributes(res) <- "MS:4000175"
    
    if (identificationLevel %in% c("all", "identified"))
        attributes(res) <- list(medianCharge = ms_term)
    
    res
}


