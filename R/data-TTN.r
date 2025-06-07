#' TTN data set
#' @name TTN
#' @aliases TTN TTN.gen TTN.fam TTN.bim TTN.pop
#' @docType data
#' 
#' @description
#' These data have been extracted from the 1000 Genomes data.
#' The data set contains the genotype matrix \code{TTN.gen}, the pedigree matrix \code{TTN.fam} and a matrix \code{TTN.bim},
#' corresponding to 503 individuals of European populations and 733 SNPs on chromosome 2, on a ~600kb segment
#' containing the Titin gene. There is also a factor \code{TTN.pop}, which gives the population from which each
#' individual is drawn (CEU = Utah residents of Northern Western European ancestry, FIN = Finnish, GBR = England and Scottland,
#' IBS = Iberian, TSI = Toscani).
#' 
#' 
#' 
#' 
#' 
#' @format
#' There are three data objects in the dataset:
#' \describe{
#' \item{list("TTN.gen")}{ Genotype matrix }
#' \item{list("TTN.fam")}{ Data frame containing all variables corresponding to a \code{.fam} file  }
#' \item{list("TTN.bim")}{ Data frame containing all variables corresponding to a \code{.bim} file }
#' \item{list("TTN.pop")}{ Factor giving the population from which each individual is drawn }
#' }
#' 
#' 
#' @references McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632
#' 
#' @source  The data were obtained from the 1000 Genomes project (see \url{https://www.internationalgenome.org/}).
#' @keywords datasets
#' @examples
#' 
#' data(TTN)
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' x
#' 
NULL



