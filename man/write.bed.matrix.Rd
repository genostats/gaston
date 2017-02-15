\name{write.bed.matrix}
\alias{write.bed.matrix}

\title{ Save a \code{\link{bed.matrix}} }

\description{
  Save a \code{\link{bed.matrix}} in several files
}
\usage{
write.bed.matrix(x, basename, bed = paste(basename, ".bed", sep=""),   
                              fam = paste(basename, ".fam", sep=""),
                              bim = paste(basename, ".bim", sep=""), 
                              rds = paste(basename, ".rds", sep=""))
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{basename}{ Basename of all files }
  \item{bed}{ Name of the \code{.bed} file }
  \item{fam}{ Name of the \code{.fam} file }
  \item{bim}{ Name of the \code{.bim} file }
  \item{rds}{ Name of the \code{.rds} file }
}

\details{ If any of \code{bed}, \code{fam}, \code{bim} and \code{rds} is \code{NULL},
the corresponding file will not be written. 

The \code{.fam} and \code{.bim} files are useful for reading files with other softwares.
The \code{.rds} file can be read by \code{read.bed.matrix}.

The \code{.bed}, \code{.fam} and \code{.bim} files follow the PLINK specifications
(\url{http://zzz.bwh.harvard.edu/plink/binary.shtml}).
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{read.bed.matrix}}, \code{\link[base:saveRDS]{saveRDS}} }

\examples{
# Load data
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

# Write object in LCT.bed and LCT.RData
\dontrun{
write.bed.matrix(x, "LCT")
}
}
