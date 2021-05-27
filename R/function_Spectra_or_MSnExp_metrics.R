
#' @name 
#' 
#' @title 
#' 
#' @description 
#' 
#' @details
#' 
#' @param
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 

#' @name rtimeDuration
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
rtimeDuration <- function(object) {
    RT <- rtime(object = )
    rtDuration <- max(RT) - min(RT)
    return(rtDuration)
}


#' @name RToverTICquantile
#' 
#' @title RT over TIC quantile (QC:4000054)
#' 
#' @description 
#' The interval when the respective quantile of the TIC accumulates divided by 
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
#' @param msnexp `MSnExp` object
#' 
#' @return numeric(4)
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
RToverTICquantile <- function(object) {
    
    if (is(object, "MSnExp"))
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
    ticSum <- cumsum(TIC)
    RT <- rtime(object)
    
    ## create relative retention time ############ assume that rt always start at 0?????????
    RT <- RT / max(RT)
    
    quantileTICSum <- quantile(ticSum)
    
    ############### my understanding -->  #############################
    ## at which observed RT event (present in object) does the 
    ## (theoretical) quantile value exceed the measured cumulative intensity, 
    ## return the maximum index
    indMax <- lapply(seq_along(quantileTICSum), 
                        function(x) max(which(ticSum <= quantileTICSum[[x]])))
    indMax <- unlist(indMax)
    
    ## return the relative retention time when quantile exceeds cumulated 
    ## (measured) intensites
    quantileRT <- RT[indMax]
    names(quantileRT) <- names(quantileTICSum)
    
    return(quantileRT)
}

#' @name MSQuantilesAlongRT
#' 
#' @title MS1/MS2 quantiles RT fraction (QC:4000055/QC:4000056)
#' 
#' @description
#' "The interval used for acquisition of the first, second, third, and fourth 
#' quarter of all MS1 events divided by RT-Duration." [PSI:QC]
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000021 ! retention time metric
#' is_a: QC:4000023 ! MS1 metric
#' 
#' @param object `MSnExp` or `Spectra` object
#' @param MSLevel `numeric(1)`
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
MSQuantilesAlongRT <- function(object, MSLevel = 1L) {
    
    ## truncate objecg based on the MSLevel
    object <- filterMsLevel(object = object, MSLevel)
    
    RT <- rtime(object) ########### is this sorted???
    ## create relative retention time ############ assume that rt always start at 0?????????
    RT <- RT / max(RT)   
    
    ## partition the object (rows) into four parts 
    ## (they are not necessarily equal)
    ind <- rep(seq_len(4), length.out = length(object))
    ind <- sort(ind)
    
    ## get the last retention time event that falls within the partition group
    rtimeGroup <- lapply(seq_len(4), function(x) max(RT[ind == x]))
    rtimeGroup <- unlist(rtimeGroup)
    
    return(rtimeGroup)
}

#' @name MSquantileTICratioToQuantiles
#' 
#' @title MS1 quantile TIC change ratio to Quantile 1 (QC:4000057) or to 
#' previous quantile (QC:4000058)
#' 
#' @description 
#' "The log ratio for the second to n-th quantile of TIC changes over first 
#' quantile of TIC changes." [PSI:QC]
#' id: QC:4000057
#' def: "The log ratio for the second to n-th quantile of TIC over the previous 
#' quantile of TIC. For the boundary elements min/max are used." [PSI:QC]
#' id: QC:4000058
#' 
#' @details
#' is_a: QC:4000004 ! n-tuple
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000022 ! chromatogram metric
#' is_a: QC:4000023 ! MS1 metric
#' 
#' The `log2` values are returned instead of the `log` values.
#' 
#' @param msnexp `MSnExp` object
#' @param relativeTo `character`
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
MSquantileTICratiotoQuantiles <- function(object, relativeTo = c("Q1", "previous")) {
    
    relativeTo <- match.arg(relativeTo)
    
    ## create cumulative sum of tic/ionCount
    if (is(object, "MSnExp"))
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
    ticSum <- cumsum(TIC)
    
    quantileTICSum <- quantile(ticSum)
    
    ## calculate the changes in TIC per quantile
    changeQ1 <- quantileTICSum[["25%"]] - quantileTICSum[["0%"]]
    changeQ2 <- quantileTICSum[["50%"]] - quantileTICSum[["25%"]]
    changeQ3 <- quantileTICSum[["75%"]] - quantileTICSum[["50%"]]
    changeQ4 <- quantileTICSum[["100%"]] - quantileTICSum[["75%"]]

    ## calculate the ratio between Q2/Q3/Q4 to Q1 TIC changes
    if (relativeTo == "Q1") 
        ratioQuantileTIC <- c(changeQ2, changeQ3, changeQ4) / changeQ1

    ## calculate the ratio between Q2/Q3/Q4 to previous quantile TIC changes
    if (relativeTo == "previous")
        ratioQuantileTIC <- c(changeQ2 / changeQ1, changeQ3 / changeQ2, 
                                                        changeQ4 / changeQ3)
    
    ## take the log2 and return
    logRatioQuantileTIC <- log2(ratioQuantileTIC)
    
    return(logRatioQuantileTIC)
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
#' @param object `MSnExp` or `Spectra` object
#' @param MSLevel `numeric(1)`
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
numberSpectra <- function(object, MSLevel = 1L) {
    
    object <- filterMsLevel(object = object, MSLevel)
    len <- length(object)
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
#' @param object `MSnexp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
medianPrecursorMZ <- function(object) {
    
    ################ FDR correction???????????
    mz <- Spectra::precursorMz(object)
    medianMZ <- median(mz, na.rm = TRUE)
    return(medianMZ)
}

#' @name rtimeIQR
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
rtimeIQR <- function(object) {
    RT <- rtime(object)
    ## IQR???, what is the unit for rtime, always seconds??
    iqr <- IQR(RT)

    return(iqr)
}

#' @name rtimeIQRrate
#' 
#' @title Peptide identification rate of the interquartile RT period (QC:4000073)
#' 
#' @description
#' The identification rate of peptides for the interquartile retention time 
#' period, in peptides per second." [PSI:QC]
#' id: QC:4000073
#' 
#' @details
#' Higher rates indicate efficient sampling and identification.
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000022 ! chromatogram metric
#' 
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(2)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
rtimeIQRrate <- function(object) {
    
    RT <- rtime(object) 
    quantileRT <- quantile(RT)
    
    ## get the RT values of the 25% and 75% quantile
    quantile25RT <- quantileRT[["25%"]]
    quantile75RT <- quantileRT[["75%"]]
    
    ## get the number of eluted features between the 25% and 75% quantile
    nFeatures <- RT >= quantile25RT & RT <= quantile75RT
    nFeatures <- sum(nFeatures)
    
    ## divide the number of eluted features between the 25% and 75% quantile
    ## by the IQR to get the elution rate per second 
    rate <- nFeatures / rtimeIQR(object)
    
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
areaUnderTIC <- function(object) {
    if (is(object, "MSnExp"))
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(4)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
areaUnderTICRTquantiles <- function(object) {
    
    if (is(object, "MSnExp"))
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
    RT <- rtime(object)
    quantileRT <- quantile(RT)
    
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
extentIdentifiedPrecursorIntensity <- function(object) {
    
    ## retrieve the precursorIntensity and calculate the 5% and 95% quantile
    precInt <- precursorIntensity(object)
    quantilePrecInt <- quantile(precInt, probs = c(0.05, 0.95), na.rm = TRUE)
    
    ## calculate the ratio between the 95% and 5% quantile and return the value
    ratio <- quantilePrecInt[["95%"]] / quantilePrecInt[["5%"]]
    
    return(ratio)
}

#' @name medianTICrtimeIQR
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
#' The function `medianTICrtimeIQR` uses the function `ionCount` as an 
#' equivalent to the TIC.
#' 
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#'
#' @author 
#'
#' @export
#'
#' @examples 
#'
medianTICrtimeIQR <- function(object, MSLevel = 1L) {
    
    object <- filterMsLevel(object = object, MSLevel)
    
    ## get the Q1 to Q3 of identifications 
    ## (half of peptides that are identitied)
    ind <- rep(seq_len(4), length.out = length(object))
    ind <- sort(ind)
    Q1ToQ3 <- object[ind %in% c(2, 3), ]

    ## take the TIC of the Q1 to Q3 of identifications
    if (is(object, "MSnExp"))
        ticQ1ToQ3 <- ProtGenerics::tic(Q1ToQ3)
    
    ## how to define TIC? take the ionCount?? #######################
    ## take the ionCount of the Q1 to Q3 of identifications
    if (is(object, "Spectra"))
        ticQ1ToQ3 <- ProtGenerics::ionCount(Q1ToQ3)
    
    ## take the median value of the TIC within this interval and return
    medianTICQ1ToQ3 <- median(ticQ1ToQ3)
    
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

#' @param object `MSnExp` or `Spectra` object
#' @param MSLevel `numeric`
#' 
#' @return
#' `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
medianTICofRTRange <- function(object, MSLevel = 1L) {
    
    object <- filterMsLevel(object = object, MSLevel)
    
    if (is(object, "MSnExp"))
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
    ## retrieve retention time
    RT <- ProtGenerics::rtime(object)
    
    ## retrieve number of features in object and calculate the number for half
    ## of the features
    n <- length(object)
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
    medianTICMin <- median(ticMin, na.rm = TRUE)
    
    return(medianTICMin)
}

#' @name MZacquisitionRange
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
#' @param object `MSnExp` or `Spectra` object
#' @param MSLevel `numeric`
#' 
#' @return
#' `numeric(2)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' MZacquisitionRange(spectra)
#' 
MZacquisitionRange <- function(object, MSLevel = 1L) {
    
    object <- filterMsLevel(object = object, MSLevel)
    
    mzList <- ProtGenerics::mz(object)
    MZ <- unlist(mzList)
    mzRange <- range(MZ)
    
    return(mzRange)
}

#' @name RTacquisitionRange
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return
#' `numeric(2)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' RTacquisitionRange(spectra)
#' 
RTacquisitionRange <- function(object) {
    RT <- ProtGenerics::rtime(object)
    rtRange <- range(RT)
    
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return
#' `numeric(2)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
#' 
precursorIntensityRange <- function(object) {
    
    int <- ProtGenerics::precursorIntensity(object)
    rangeInt <- range(int)
    
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(3)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples
#' precursorIntensityQuartiles(object)
precursorIntensityQuartiles <- function(object) {
    
    int <- ProtGenerics::precursorIntensity(object)
    quartiles <- quantile(int, probs = c(0.25, 0.50, 0.75))
    
    return(quartiles)
}


#' @name precursorIntensityMean
#' 
#' @title Precursor intensity distribution mean (QC:4000168),
#' Identified precursor intensity distribution mean (QC:4000229), or
#'  Unidentified precursor intensity distribution mean (QC:4000234)
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples
#' precursorIntensityMean(object)
precursorIntensityMean <- function(object) {
    int <- ProtGenerics::precursorIntensity(object)
    intMean <- mean(int, na.rm = TRUE)
    
    return(mzMean)
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
#' @param object `MSnExp` or `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples
#' precursorIntensitySD(object)
precursorIntensitySD <- function(object) {
    int <- ProtGenerics::precursorIntensity(object)
    mzSD <- sd(int, na.rm = TRUE)
    
    return(mzSD)
}

#' @name MS1Signal10XChange
#' 
#' @title MS1 signal jump/fall (10x) count (QC:4000172/QC:400173)
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
#' @param object `MSnExp` or `Spectra` object
#' @param change `character`, one of `"jump"` or `"fall"`
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
MS1Signal10XChange <- function(object, change = c("jump", "fall")) {
    
    change <- match.arg(change)
    
    ########## does this make sense only for MSnExp?????
    if (is(object, "MSnExp")) 
        TIC <- tic(object)
    if (is(object, "Spectra"))
        TIC <- ionCount(object)
    
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

