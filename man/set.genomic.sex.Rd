\name{set.genomic.sex}
\alias{set.genomic.sex}

\title{ Genomic Sex }

\description{
  Returns an updated \code{\link{bed.matrix}} with a new variable for the genomic sex
  of each individual.
}
\usage{ set.genomic.sex(x, plot = FALSE, verbose = getOption("gaston.verbose",TRUE)) }
\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{plot}{ If \code{TRUE}, plots the variables used for the clustering }
  \item{verbose}{ If \code{TRUE}, displays information on the function actions }
}

\details{
  For each individual, the function uses the hetorozygosity rate for SNPs on X chromosome, 
  and the call rate for SNPs on the Y chromosomes (both statistics computed by \code{\link{set.stats}}), 
  to cluster the individuals using \code{\link{kmeans}}. 

  If \code{plot = TRUE}, a plot is produced with the two variables used and the clusters
  determined by \code{\link{kmeans}}.
}

\value{
  A \code{\link{bed.matrix}} similar to \code{x}, with a new variable \code{x@ped$genomic.sex}
  containing the genomic sex for each individual.
}


\author{Hervé Perdry}

\seealso{ \code{\link{set.stats}}, \code{\link{set.hwe}} }

\keyword{ Genomic Sex }

