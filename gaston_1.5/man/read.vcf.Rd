\name{read.vcf}
\alias{read.vcf}

\title{ Create a \code{\link{bed.matrix}} from VCF files }

\description{ Create a \code{\link{bed.matrix}} from a \code{.vcf} file.}
\usage{ read.vcf(filename, max.snps, verbose = getOption("gaston.verbose",TRUE)) }

\arguments{
  \item{filename}{ Name of VCF file }
  \item{max.snps}{ The maximal number of SNPs to read. }
  \item{verbose}{ If \code{TRUE}, display information on the function progress }
}

\details{
This function relies on the package \code{\link[WhopGenome:WhopGenome-package]{WhopGenome}}.
Details about vcf format are given in \url{http://www.1000genomes.org/wiki/Analysis/Variant\%20Call\%20Format/vcf-variant-call-format-version-41}
}

\value{ A \code{\link{bed.matrix}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{read.bed.matrix}} }

\examples{
# Read vcf file
filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
x <- read.vcf( filepath )
x
}
\keyword{ vcf files }
