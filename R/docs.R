##' A package to perform robust quantitative traits association testing of copy number variants. It provides two methods for association study: first, the observed probe intensity measurement can be directly used to detect the association of CNV with phenotype of interest. Second, the most probable copy number is estimated with the proposed likelihood and the association of the most probable copy number with phenotype is tested. Also, it can be used to determine the optimal clustering number and clustering assignment for each individuals. This method can be applied to both the independent and correlated population.
##'
##' \tabular{ll}{
##' Package: \tab PedCNV\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1\cr
##' Date: \tab 2013-09-03\cr
##' License: \tab MIT \cr
##' Main functions:
##' \tab AssoTestProc \cr
##' \tab ClusProc \cr
##' \tab STE \cr
##' \tab STIM \cr
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

##' Simulated data. It contains CNV intensity measurements, environment variables and FAM file which follows the same format defined in PLINK.
##' 
##' @title CNV simulated data
##' @name simudat
##' @docType data
##' @author Meiling Liu
##' @examples
##' data(simudat)
NULL

##' Empirical/kinship correlation matrix between individuals. This correlation matrix can be calculated based on the familial relationship between individuals or large-scale SNP data by omic data analysis toolkit WISARD. The free software WISARD can be downloaded from \url{http://biostat.cau.ac.kr/wisard/}. If correlation matrix is estimated with the large-scale SNP data, the proposed method becomes robust under the presence of population substructure.
##' 
##' @name phi
##' @title Empirical correlation matrix
##' @docType data
##' @author Meiling Liu
##' @references WISARD \url{http://biostat.cau.ac.kr/wisard/}
##' @examples
##' data(phi)
NULL


