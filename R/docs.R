##' A package to perform robust quantitative traits association testing of copy number variants.
##' 
##' @name PedCNV-package
##' @docType package
##' @title PedCNV: CNV association studies
##' @author Meiling \email{meiling.sta@@gmail.com}
##' @keywords Association; Clustering; CNV
##' @useDynLib PedCNV
##' @references On the association analysis of CNV data: fast and efficient method with family-based samples
NULL

##' Simulated data set. It contains CNV intensity measurements, environment variables and PED file which follows the same format defined in PLINK.
##' 
##' @title CNV simulated data
##' @name dat
##' @format A list.
##' @docType data
##' @author Meiling \email{meiling.sta@@gmail.com}
##' @usage data(dat)
NULL

##' Empirical correlation matrix. This correlation matrix is calculated based on the family structure in dat by using omic data analysis toolkit WISARD. This toolkit will eventually included in PedCNV.
##' 
##' @name phi
##' @title Empirical correlation matrix
##' @format A matrix.
##' @docType data
##' @author Meiling \email{meiling.sta@@gmail.com}
##' @references WISARD \url{http://biostat.cau.ac.kr/wisard/}
##' @usage data(phi)
NULL


