#' @name Lee_2019_meta_vals
#' 
#' @aliases meta
#' @aliases vals
#' 
#' @title Example data for \code{MsQuality}: data set of Lee et al. (2019)
#' 
#' @description 
#' The data set of Lee et al. (2019) contains metabolite information measured
#' by reverse phase liquid chromatography (RPLC) coupled to mass spectrometry 
#' and hydrophilic interaction liquid chromatography (HILIC) coupled to mass 
#' spectrometry (file `STables - rev1.xlsx` in the Supplementary Information).
#' The xlsx sheets `Methods` and `Raw data` were stored as txt files.
#' 
#' \code{Lee_2019_meta_vals} contains two data frame objects:
#' one containing information on metabolite meta-data and one containing 
#' intensity values on metabolites. 
#' The object will be used as an example data set in the vignette to
#' show the functionality of the packages.
#' 
#' @references 
#' Lee et al. (2019). A large-scale analysis of targeted metabolomics data from 
#' heterogeneous biological samples provides insights into metabolite dynamics.
#' Metabolomics, 103, doi: 10.1007/s11306-019-1564-8. 
#' 
#' @docType data
#' 
#' @return \code{data.frame}
#'
#' @format \code{data.frame}
#'
#' @source
#' path_to_meta <- "Lee_et_al_2019_Stables_rev1_Methods.txt"
#' meta <- read.delim(path_to_meta, dec = ".", header = TRUE)
#' 
#' ## print number of metabolites per measurement (meta data)
#' table(meta$Method)
#' 
#' path_to_vals <- "Lee_et_al_2019_Stables_rev1_Raw_data.txt"
#' vals <- read.delim(path_to_vals, dec = ".", header = TRUE)
#' 
#' ## print number of metabolites per measurement (intensity data)
#' table(grepl(vals$Metabolite, pattern = "_rp$"))
#' table(grepl(vals$Metabolite, pattern = "_hn$"))
#' 
#' ## save the two objects as an RData object
#' save(meta, vals, file = "Lee_2019_meta_vals.RData", compress = "xz")
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
