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


library(Spectra)

spd <- DataFrame(
    msLevel = c(2L, 2L, 2L, 2L, 2L, 2L),
    polarity = c(1L, 1L, 1L, 1L, 1L, 1L),
    id = c("HMDB0000001", "HMDB0000001", "HMDB0001847", "HMDB0000517)", "HMDB0000177", "HMDB0000162"),
    name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine", "Arginine", "Histidine", "Proline"))

## Assign m/z and intensity values.
spd$mz <- list(
    c(109.2, 124.2, 124.5, 170.16, 170.52),
    c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
    c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
      111.0551, 123.0429, 138.0662, 195.0876),
    c(60.00, 70.00, 116),
    c(109.12, 156.07),
    c(69.91, 116.07))
spd$intensity <- list(
    c(3.407, 47.494, 3.094, 100.0, 13.240),
    c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
    c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994),
    c(200, 1000, 290),
    c(300, 1100),
    c(202, 471))
sps <- Spectra(spd)
sps$rtime <- c(1, 3, 7, 9, 10, 13)
sps
spectra <- sps

library(xcms)
swath_file <- system.file("TripleTOF-SWATH",
                          "PestMix1_SWATH.mzML",
                          package = "msdata")
swath_data <- readMSData(swath_file, mode = "onDisk")
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaks(swath_data, param = cwp)
cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp)
swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.9,
                                             return.type = "Spectra")


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
#' @param
#' 
#' @return numeric(4)
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
RToverTICquantile <- function(spectra) {
    spectraIonCount <- ionCount(spectra) ## or precursorIntensity???????????
    spectraIonCountSum <- cumsum(spectraIonCount)
    spectraRtime <- rtime(spectra)
    
    ## create relative retention time ############ assume that rt always start at 0?????????
    spectraRtime <- spectraRtime / max(spectraRtime)
    
    quantileIonCountSum <- quantile(spectraIonCountSum)
    
    ############### my understanding -->  #############################
    ## at which observed RT event (present in spectra) does the 
    ## (theoretical) quantile value exceed the measured cumulative intensity, 
    ## return the maximum index
    indMax <- lapply(seq_along(quantileIonCountSum), 
        function(x) max(which(spectraIonCountSum <= quantileIonCountSum[[x]])))
    indMax <- unlist(indMax)
    
    ## return the relative retention time when quantile exceeds cumulated 
    ## (measured) intensites
    quantileRtime <- spectraRtime[indMax]
    names(quantileRtime) <- names(quantileIonCountSum)
    
    return(quantileRtime)
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
#' @param spectra `SpectraObject`
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
MSQuantilesAlongRT <- function(spectra, MSLevel = 1) {
 
    ## truncate spectra based on the MSLevel
    spectra <- spectra[msLevel(spectra) == MSLevel, ]
    
    spectraRtime <- rtime(spectra)
    ## create relative retention time ############ assume that rt always start at 0?????????
    spectraRtime <- spectraRtime / max(spectraRtime)   
    
    ## partition the spectra (rows) into four parts 
    ## (they are not necessarily equal)
    ind <- rep(seq_len(4), length.out = length(spectra))
    ind <- sort(ind)
    
    ## get the last retention time event that falls within the partition group
    rtimeGroup <- lapply(seq_len(4), function(x) max(spectraRtime[ind == x]))
    rtimeGroup <- unlist(rtimeGroup)
    
    return(rtimeGroup)
}

#' @name MSquantileTICratioToQuantiles
#' 
#' @title MS1 quantile TIC change ratio to Quantile 1 (QC:4000057) or to 
#' previous quantile (QC:4000058)
#' 
#' @description 
#' "The log ratio for the second to n-th quantile of TIC changes over first quantile of TIC changes." [PSI:QC]
#' id: QC:4000057
#' def: "The log ratio for the second to n-th quantile of TIC over the previous quantile of TIC. For the boundary elements min/max are used." [PSI:QC]
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
#' @param spectra `Spectra` object
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
MSquantileTICratiotoQuantiles <- function(spectra, relativeTo = c("Q1", "previous")) {
    
    relativeTo <- match.arg(relativeTo)
    
    spectraIonCount <- ionCount(spectra) ## or precursorIntensity???????????
    spectraIonCountSum <- cumsum(spectraIonCount)
    
    quantileIonCountSum <- quantile(spectraIonCountSum)
    
    ## calculate the changes in TIC per quantile
    changeQ1 <- quantileIonCountSum[["25%"]] - quantileIonCountSum[["0%"]]
    changeQ2 <- quantileIonCountSum[["50%"]] - quantileIonCountSum[["25%"]]
    changeQ3 <- quantileIonCountSum[["75%"]] - quantileIonCountSum[["50%"]]
    changeQ4 <- quantileIonCountSum[["100%"]] - quantileIonCountSum[["75%"]]
    
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
#' @param spectra `Spectra` object
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
numberSpectra <- function(spectra, MSLevel = 1) {
    
    len <- sum(msLevel(spectra) == MSLevel)
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
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
medianPrecursorMZ <- function(spectra) {
    mz <- Spectra::precursorMz(spectra)
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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
rtimeIQR <- function(spectra) {
    spectraRT <- rtime(spectra)
    ## IQR???, what is the unit for rtime, always seconds??
    iqr <- IQR(spectraRT)
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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(2)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
rtimeIQRrate <- function(spectra) {
    spectraRT <- rtime(spectra)
    quantileSpectraRT <- quantile(spectraRT)
    
    ## get the RT values of the 25% and 75% quantile
    quantile25RT <- quantileSpectraRT[["25%"]]
    quantile75RT <- quantileSpectraRT[["75%"]]
    
    ## get the number of eluted features between the 25% and 75% quantile
    nFeatures <- spectraRT >= quantile25RT & spectraRT <= quantile75RT
    nFeatures <- sum(nFeatures)
    
    ## divide the number of eluted features between the 25% and 75% quantile
    ## by the IQR to get the elution rate per second 
    rate <- nFeatures / rtimeIQR(spectra)
    
    return(rate)
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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(1)`
#'
#' @author 
#'
#' @export
#'
#' @examples 
#'
medianTICrtimeIQR <- function(spectra) {
    
    ## get the Q1 to Q3 of identifications 
    ## (half of peptides that are identitied)
    ind <- rep(seq_len(4), length.out = length(spectra))
    ind <- sort(ind)
    spectraQ1ToQ3 <- spectra[ind %in% c(2, 3), ]
    
    ## how to define TIC? take the ionCount?? #######################
    ## take the ionCount of the Q1 to Q3 of identifications and calculate the
    ## TIC by adding up the individual intensities per intensity
    ionCountQ1ToQ3 <- ionCount(spectraQ1ToQ3)
    
    ## take the median value of the ionCount within this interval and return
    medianIonCountQ1ToQ3 <- median(ionCountQ1ToQ3)
    
    return(ionCountQ1ToQ3)
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

#' @param spectra `Spectra` object
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
MZacquisitionRange <- function(spectra) {
    spectraMZList <- mz(spectra)
    spectraMZ <- unlist(spectraMZList)
    range(spectraMZ)
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
#' @param spectra `Spectra` object
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
RTacquisitionRange <- function(spectra) {
    spectraRT <- rtime(spectra)
    range(spectraRT)
}

#' @name peakDensityMS1Mean
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


#' @name peakDensityMS1SD
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


#' @name precursorIntensityMean
#' 
#' @title Precursor intensity distribution mean (QC:4000168)
#' 
#' @description 
#' "From the distribution of precursor intensities, the mean." [PSI:QC]
#' id: QC:4000168
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples
#' precursorIntensityMean(spectra)
precursorIntensityMean <- function(spectra) {
    mz <- precursorMz(spectra)
    mean(mz)
}

#' @name precursorIntensitySD
#' 
#' @title Precursor intensity distribution sigma (QC:4000169)
#' 
#' @description 
#' "From the distribution of precursor intensities, the sigma value." [PSI:QC]
#' id: QC:4000169
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000010 ! ID free
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra
#' 
#' @return `numeric(1)`
#' 
#' @author 
#' 
#' @export
#' 
#' @examples
#' precursorIntensitySD(spectra)
precursorIntensitySD <- function(spectra) {
    mz <- precursorMz(spectra)
    sd(mz)
}

#' @name RatioCharge1over2
#' 
#' @title Charged peptides ratio 1+ over 2+ (QC:4000174)
#' 
#' @description 
#' "Ratio of 1+ peptide count over 2+ peptide count in identified spectra" [PSI:QC]
#' id: QC:4000174
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` objects
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
RatioCharge1over2 <- function(spectra) {
    ############# is count == number or intensity? #############################
    charge <- Spectra::precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(1, 2) %in% chargeTable))
        chargeRatio <- chargeTable[["1"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name RatioCharge3over2
#' 
#' @title Charged peptides ratio 3+ over 2+ (QC:4000175)
#' 
#' @description 
#' "Ratio of 3+ peptide count over 2+ peptide count in identified spectra" [PSI:QC]
#' id: QC:4000175
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` objects
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
RatioCharge3over2 <- function(spectra) {
    ############# is count == number or intensity? #############################
    charge <- Spectra::precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(2, 3) %in% chargeTable))
        chargeRatio <- chargeTable[["3"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}

#' @name RatioCharge4over2
#' 
#' @title Charged peptides ratio 4+ over 2+ (QC:4000176)
#' 
#' @description 
#' "Ratio of 4+ peptide count  over 2+ peptide count  in identified spectra" [PSI:QC]
#' id: QC:4000176
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` objects
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
RatioCharge4over2 <- function(spectra) {
    ############# is count == number or intensity? #############################
    charge <- Spectra::precursorCharge(spectra)
    
    ## get the number of precursor per charge
    chargeTable <- table(charge)
    
    if (all(c(2, 4) %in% chargeTable))
        chargeRatio <- chargeTable[["4"]] / chargeTable[["2"]]
    else 
        chargeRatio <- NA
    
    return(chargeRatio)
}


#' @name meanCharge
#' 
#' @title Mean charge in identified spectra (QC:4000177)
#' 
#' @description 
#' "Mean charge in identified spectra" [PSI:QC]
#' id: QC:4000177
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
meanCharge <- function(spectra) {
    charge <- Spectra::precursorCharge(spectra)
    chargeMean <- mean(charge, na.rm = TRUE)
    return(chargeMean)
}

#' @name medianCharge
#' 
#' @title Median charge in identified spectra (QC:4000178)
#' 
#' @description 
#' "Median charge in identified spectra" [PSI:QC]
#' id: QC:4000178
#' 
#' @details
#' is_a: QC:4000003 ! single value
#' is_a: QC:4000009 ! ID based
#' is_a: QC:4000001 ! QC metric
#' 
#' @param spectra `Spectra` object
#' 
#' @return
#' 
#' @author 
#' 
#' @export
#' 
#' @examples 
#' 
medianCharge <- function(spectra) {
    charge <- Spectra::precursorCharge(spectra)
    chargeMedian <- median(charge, na.rm = TRUE)
    return(chargeMedian)
}


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