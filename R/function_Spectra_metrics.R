
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


#' @name RatioCharge1over2
#' 
#' @title Charged peptides ratio 1+ over 2+ (QC:4000174) or 
#' Charged spectra ratio +1 over +2 (QC:4000179)
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
    ## is there a way to get charge of actual entries, not only of precursor?
    ############# is count == number or intensity? #############################
    charge <- ProtGenerics::precursorCharge(spectra)
    
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
#' @title Charged peptides ratio 3+ over 2+ (QC:4000175) or 
#' charged spectra ratio +3 over +2 (QC:4000180)
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
    ## is there a way to get charge of actual entries, not only of precursor?
    ############# is count == number or intensity? #############################
    charge <- ProtGenerics::precursorCharge(spectra)
    
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
    ## is there a way to get charge of actual entries, not only of precursor?
    ############# is count == number or intensity? #############################
    charge <- ProtGenerics::precursorCharge(spectra)
    
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
    charge <- ProtGenerics::precursorCharge(spectra)
    chargeMedian <- median(charge, na.rm = TRUE)
    return(chargeMedian)
}

