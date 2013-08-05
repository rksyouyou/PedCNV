##' A package to perform robust quantitative traits association testing of copy number variants. It provides two methods for association study: first, the observed probe intensity measurement can be directly used to detect the association of CNV with phenotype of interest. Second, the most probable copy number is estimated with the proposed likelihood and the association of the most probable copy number with phenotype is tested. Also, it can be used to determine the optimal clustering number and clustering assignment for each individuals. This method can be applied to both the independent and correlated population.
##'
##' \tabular{ll}{
##' Package: \tab PedCNV\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1\cr
##' Date: \tab 2013-08-03\cr
##' License: \tab MIT \cr
##' Main functions:
##' \tab AssoTestProc \cr
##' \tab ClusProc \cr
##' \tab print.asso\cr
##' \tab print.clus\cr
##' \tab plot.clus\cr
##' }
##'
##' @name PedCNV-package
##' @docType package
##' @title CNV association implementation
##' @author Meiling Liu, Sungho Won and Weicheng Zhu
##' @useDynLib PedCNV
##' @references On the association analysis of CNV data: fast and efficient method with family-based samples
NULL

##' Simulated data. It contains CNV intensity measurements, environment variables and PED file which follows the same format defined in PLINK.
##' 
##' @title CNV simulated data
##' @name simudat
##' @docType data
##' @author Meiling Liu
##' @examples
##' data(simudat)
NULL

##' Empirical correlation matrix. This correlation matrix is calculated based on the family structure in simudat by using omic data analysis toolkit WISARD. This toolkit will eventually be included in PedCNV.
##' 
##' @name phi
##' @title Empirical correlation matrix
##' @docType data
##' @author Meiling 
##' @references WISARD \url{http://biostat.cau.ac.kr/wisard/}
##' @examples
##' data(phi)
NULL


