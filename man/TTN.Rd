\name{TTN}
\docType{data}
\alias{TTN}
\alias{TTN.gen}
\alias{TTN.fam}
\alias{TTN.bim}
\alias{TTN.pop}

\title{TTN data set}

\description{
  These data have been extracted from the 1000 Genomes data.
  The data set contains the genotype matrix \code{TTN.gen}, the pedigree matrix \code{TTN.fam} and a matrix \code{TTN.bim},
  corresponding to 503 individuals of European populations and 733 SNPs on chromosome 2, on a ~600kb segment 
  containing the Titin gene. There is also a factor \code{TTN.pop}, which gives the population from which each
  individual is drawn (CEU = Utah residents of Northern Western European ancestry, FIN = Finnish, GBR = England and Scottland, 
  IBS = Iberian, TSI = Toscani).
}
\usage{ data(TTN) }

\format{
  There are three data objects in the dataset:
  \describe{
    \item{\code{TTN.gen}}{ Genotype matrix }
    \item{\code{TTN.fam}}{ Data frame containing all variables corresponding to a \code{.fam} file  }
    \item{\code{TTN.bim}}{ Data frame containing all variables corresponding to a \code{.bim} file }
    \item{\code{TTN.pop}}{ Factor giving the population from which each individual is drawn }
  }
}

\source{ The data were obtained from the 1000 Genomes project (see \url{https://www.internationalgenome.org/}). }

\references{McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }

\examples{
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
x
}

\keyword{datasets}
