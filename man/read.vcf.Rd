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
## Read vcf file (from file name)
filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
x1 <- read.vcf( filepath )
x1

## from file handle given by WhopGenome 
# open file
handle <- WhopGenome::vcf_open(filepath)
# list all samples
samples <- WhopGenome::vcf_getsamples(handle)
# keep only samples 1, 7 and 10
WhopGenome::vcf_selectsamples(handle, samples[c(1,7,10)] )
# keep only snps on chr 2 from pos = 136400000 to pos = 136500000
WhopGenome::vcf_setregion(handle, "2", 136400000, 136500000)
# make bed matrix from data
x2 <- read.vcf(handle)
# and close the file
WhopGenome::vcf_close(handle)
x2

# we can check that this is the same as extracting a part of x1
all( as.matrix(x1[c(1,7,10), 1:207]) == as.matrix(x2) )

}
\keyword{ vcf files }
