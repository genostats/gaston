\name{AGT}
\docType{data}
\alias{AGT}
\alias{AGT.gen}
\alias{AGT.fam}
\alias{AGT.bim}
\alias{AGT.pop}

\title{AGT data set}

\description{
  These data have been extracted from the 1000 Genomes data.
  The data set contains the genotype matrix \code{AGT.gen}, the pedigree matrix \code{AGT.fam} and a matrix \code{AGT.bim},
  corresponding to 503 individuals of European populations and 361 SNPs on chromosome 1, on a ~100kb segment 
  containing the Angiotensinogen gene. There is also a factor \code{AGT.pop}, which gives the population from which each
  individual is drawn (CEU = Utah residents of Northern Western European ancestry, FIN = Finnish, GBR = England and Scottland, 
  IBS = Iberian, TSI = Toscani).
}
\usage{ data(AGT) }

\format{
  There are three data objects in the dataset:
  \describe{
    \item{\code{AGT.gen}}{ Genotype matrix }
    \item{\code{AGT.fam}}{ Data frame containing all variables corresponding to a \code{.fam} file  }
    \item{\code{AGT.bim}}{ Data frame containing all variables corresponding to a \code{.bim} file }
    \item{\code{AGT.pop}}{ Factor giving the population from which each individual is drawn }
  }
}

\source{ The data were obtained from the 1000 Genomes project (see \url{https://www.internationalgenome.org/}). }

\references{McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }

\examples{
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
x
}

\keyword{datasets}
