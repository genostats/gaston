#' Reshape a Genetic Relationship Matrix
#' 
#' @description
#' Reshapes a GRM into a data frame listing relationship of (possibly all) pairs of individuals.
#' Options are provided to specify ranges of relationship values to include or exclude. This
#' is useful in the Quality Control process.
#' 
#' @param K  A symmetric matrix (such as produced by \code{\link{GRM}})
#' @param include  Range of values to include (default is to include all values)
#' @param exclude  Range of values to exclude (default it to exclude nothing)
#' 
#' @details
#' The relationship between individuals \eqn{i} and \eqn{j} is the coefficient \eqn{k_{ij}}{k_ij}
#' in the matrix \eqn{K}. The functions lists all pair \eqn{i, j} with \eqn{i < j} and \eqn{k_{ij}}{k_ij}
#' in the range defined by \code{include} and outside the range defined by \code{exclude}.
#' 
#' @return
#' A data frame with three columns named \code{i}, \code{j}, \code{k}.
#' @seealso  \code{\link{GRM}}
#' 
#' @keywords  Genetic Relationship Matrix
#' @examples
#' 
#' # load chr2 data set (~10k SNPs in low LD)
#' x <- read.bed.matrix( system.file("extdata", "chr2.bed", package="gaston") )
#' 
#' # Compute Genetic Relationship Matrix
#' K <- GRM(x)
#' 
#' # List all pairs if individuals with a relationship above 0.07
#' pairs <- reshape.GRM(K, exclude = c(-Inf, 0.07))
#' 
#' # Exclude first individual from each such pair
#' x1 <- x[ -pairs$i, ]
#' 
#' @export reshape.GRM
reshape.GRM <- function(K, include = c(-Inf, +Inf), exclude) {
  diag(K) <- NA
  if(missing(exclude))
    w <- which(include[1] < K & K < include[2])
  else 
    w <- which(include[1] < K & K < include[2] & (K < exclude[1] | K > exclude[2]))
  I <- row(K)[w]
  J <- col(K)[w]
  R <- K[w]
  ww <- (I < J)
  i <- I[ww];
  j <- J[ww];
  data.frame(i = i, j = j, id_i = rownames(K)[i], id_j = colnames(K)[j], k = R[ww])
}

