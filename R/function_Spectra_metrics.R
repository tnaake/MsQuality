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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge
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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge
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
#' @param spectra `Spectra` object
#' 
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge
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
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge
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
#' @return `numeric(1)`
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' 
#' @export
#' 
#' @importFrom ProtGenerics precursorCharge
#' 
#' @examples 
#' 
medianCharge <- function(spectra) {
    charge <- ProtGenerics::precursorCharge(spectra)
    chargeMedian <- median(charge, na.rm = TRUE)
    return(chargeMedian)
}

