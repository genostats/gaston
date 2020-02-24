\name{read.vcf}
\alias{read.vcf}

\title{ Create a \code{\link{bed.matrix}} from VCF files }

\description{ Create a \code{\link{bed.matrix}} from a \code{.vcf} file.}
\usage{ read.vcf(file, max.snps, get.info = FALSE, convert.chr = TRUE, 
                 verbose = getOption("gaston.verbose",TRUE)) }

\arguments{
  \item{file}{ The name of the VCF file to read}
  \item{max.snps}{ The maximal number of SNPs to read }
  \item{get.info}{ If \code{TRUE}, the INFO field from the VCF file will integrated
                   in \code{@ped$info}}
  \item{convert.chr}{ If \code{TRUE}, chromosomes ids \code{"X"}, \code{"Y"} and \code{"MT"}
                      will be converted in their numeric equivalents}
  \item{verbose}{ If \code{TRUE}, display information on the function progress }
}

\details{
The vcf format is described in \url{https://github.com/samtools/hts-specs}

In addition to the usual data in the slot \code{@snps}, the bed.matrices produced by \code{read.vcf} have
\code{@snps$quality} and \code{@snps$filter} columns corresponding to the QUAL and FILTER fields in the VCF
file. If \code{get.info = TRUE}, an additionnal column \code{@snps$info} is added, corresponding to the
INFO field.

The information about individuals in VCF files is incomplete: in the slot \code{@ped}, the columns
\code{@ped$famid} and \code{@ped$id} will both contain the sample id; sex and phenotypes will be set
to unknown.

The function currently assumes that the \code{GT} field is the first field in the genotypes format.
If it is not the case, the variants are discarded.
}

\value{ A \code{\link{bed.matrix}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{read.bed.matrix}} }

\examples{
## Read vcf file (from file name)
filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
x1 <- read.vcf( filepath )
x1
}
\keyword{ vcf files }
