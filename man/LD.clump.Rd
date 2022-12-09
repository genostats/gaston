\name{LD.clump}
\alias{LD.clump}

\title{ LD clumping }

\description{ Construct group of SNPs in LD with 'top associated SNPs' }
\usage{
LD.clump(x, p, r2.threshold, p.threshold, max.dist = 500e3)
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{p}{ A vector of p-values, or a data frame including p-values, such as sent back by \code{\link{association.test}} } 
  \item{r2.threshold}{ The maximum LD (measured by \eqn{r^2}) between SNPs in a group }
  \item{p.threshold}{ The threshold used to define associated SNPs }
  \item{max.dist}{ The maximum distance for which the LD is computed }
}

\details{
  The p-values provided through argument \code{p} are assumed to correspond to the result of an association test with the SNPs of \code{x}.

  The aim of the function is to construct cluster of SNPs in strong LD with associated SNPs.
  The algorithm first seeks the SNP with the lowest p-value (below \code{p.threshold}) ; this SNP will be the 'index' of a cluster. 
  The corresponding cluster is constructed by aggregating SNPs that are in LD (above \code{r2.threshold}) with the index. The cluster's name 
  is the position of the index SNP.
  The processus is repeated on the SNPs which are not yet attributed to a cluster, until there is no associated SNP 
  (ie SNP with a p-value below \code{threshold}) left.
  The remaining SNPs are attributed to cluster 0. 

  The LD is computed only for SNP pairs for which distance is inferior to \code{max.dist}, expressed in number of bases: above this 
  distance it is assumed to be null.
}

\value{
  If \code{p} was a data frame, then the function returns the same data frame with to extra columns, \code{cluster} and \code{is.index}.
  If \code{p} was a vector of p-values, it returns a data frame with columns \code{chr}, \code{id}, \code{pos}, \code{p}, \code{cluster} 
  and \code{is.index}.
}

\seealso{ \code{\link{LD}}, \code{\link{LD.thin}} }

\examples{
# Construct a bed matrix
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
standardize(x) <- "p"
     
# simulate quantitative phenotype with effect of SNPs #108 and #631
beta <- numeric(ncol(x))
beta[c(108,631)] <- 0.5
set.seed(1)
y <- x \%*\% beta + rnorm(nrow(x))
     
# association test with linear model 
test <- association.test(x, y, method="lm", response = "quanti")

# LD clumping
test <- LD.clump(x, test, r2.threshold = 0.25, p.threshold = 1e-8)

# use as.factor for a quick-and-dirty cluster colouring on the manhattan plot 
manhattan(test, col = as.factor(test$cluster), pch = 20)
}

\keyword{ Linkage Disequilibrium }
