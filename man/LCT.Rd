\name{LCT}
\docType{data}
\alias{LCT}
\alias{LCT.gen}
\alias{LCT.fam}
\alias{LCT.bim}
\alias{LCT.pop}

\title{LCT data set}

\description{
  These data have been extracted from the 1000 Genomes data.
  The data set contains the genotype matrix \code{LCT.gen}, the pedigree matrix \code{LCT.fam} and a matrix \code{LCT.bim},
  corresponding to 503 individuals of European populations and 607 SNPs on chromosome 2, on a ~300kb segment 
  containing the Lactase gene. There is also a factor \code{LCT.pop}, which gives the population from which each
  individual is drawn (CEU = Utah residents of Northern Western European ancestry, FIN = Finnish, GBR = England and Scottland, 
  IBS = Iberian, TSI = Toscani).

  Note that the SNP rs4988235 is associated with lactase persistence / lactose intolerence.
}
\usage{ data(LCT) }

\format{
  There are three data objects in the dataset:
  \describe{
    \item{\code{LCT.gen}}{ Genotype matrix }
    \item{\code{LCT.fam}}{ Data frame containing all variables corresponding to a \code{.fam} file  }
    \item{\code{LCT.bim}}{ Data frame containing all variables corresponding to a \code{.bim} file }
    \item{\code{LCT.pop}}{ Factor giving the population from which each individual is drawn }
  }
}

\source{ The data were obtained from the 1000 Genomes project (see \url{https://www.internationalgenome.org/}). }

\references{McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }

\examples{
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
x
which(x@snps$id == "rs4988235")
}

\keyword{datasets}
