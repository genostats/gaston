\name{read.vcf}
\alias{read.vcf}

\title{ Create a \code{\link{bed.matrix}} from VCF files }

\description{ Create a \code{\link{bed.matrix}} from a \code{.vcf} file.}
\usage{ read.vcf(file, max.snps, get.info = FALSE, convert.chr = TRUE, 
                 verbose = getOption("gaston.verbose",TRUE)) }

\arguments{
  \item{file}{ Path of VCF file or \code{WhopGenome} VCF file handle }
  \item{max.snps}{ The maximal number of SNPs to read. }
  \item{get.info}{ If \code{TRUE}, the INFO field from the VCF file will integrated
                   in \code{@ped$info}}
  \item{convert.chr}{ If \code{TRUE}, chromosomes ids \code{"X"}, \code{"Y"} and \code{"MT"}
                      will be converted in their numeric equivalents}
  \item{verbose}{ If \code{TRUE}, display information on the function progress }
}

\details{
This function relies on the package \code{\link[WhopGenome:WhopGenome-package]{WhopGenome}}.
The vcf format is described in \url{http://www.1000genomes.org/wiki/Analysis/Variant\%20Call\%20Format/vcf-variant-call-format-version-41}

In addition to the usual data in the slot \code{@snps}, the bed.matrices produced by \code{read.vcf} have
\code{@snps$quality} and \code{@snps$filter} columns corresponding to the QUAL and FILTER fields in the VCF
file. If \code{get.info = TRUE}, an additionnal column \code{@snps$info} is added, corresponding to the
INFO field.

The information about individuals in VCF files is incomplete: in the slot \code{@ped}, the columns
\code{@ped$famid} and \code{@ped$id} will both contain the sample id; sex and phenotypes will be set
to unknown.

If \code{file} is the path to a VCF file, the function will read all the samples and all the di-allelic
variants. The filtering rules provided in the \code{WhopGenome} package allows can be used by 
letting \code{file} be a filehandle obtained from \code{WhopGenome::open_vcf} (see the Examples below).

The function currently assumes that the \code{GT} field is the first field in the genotypes format.
If it is not the case, the variants are discarded.
}

\value{ A \code{\link{bed.matrix}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{read.bed.matrix}}, \code{\link[WhopGenome]{WhopGenome-package}} }

\examples{
## Read vcf file (from file name)
filepath <-system.file("extdata", "LCT.vcf.gz", package="gaston")
x1 <- read.vcf( filepath )
x1

## Read vcf from file handle given by WhopGenome 
handle <- WhopGenome::vcf_open(filepath)      # 1) open file
samples <- WhopGenome::vcf_getsamples(handle) # 2) get list of samples
# 3) select samples (here samples 1, 7 and 10)
WhopGenome::vcf_selectsamples(handle, samples[c(1,7,10)] )
# 4) [optional] select snps (here, chr 2 from pos = 136400000 to pos = 136500000)
WhopGenome::vcf_setregion(handle, "2", 136400000, 136500000)
x2 <- read.vcf(handle)                        # 5) call read.vcf
WhopGenome::vcf_close(handle)                 # 6) close file
x2

# we can check that this is the same as extracting a part of x1
all( as.matrix(x1[c(1,7,10), 1:207]) == as.matrix(x2) )

}
\keyword{ vcf files }
