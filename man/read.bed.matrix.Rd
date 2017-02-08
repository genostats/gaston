\name{read.bed.matrix}
\alias{read.bed.matrix}

\title{ Read a \code{\link{bed.matrix}} }

\description{
  Create a \code{\link{bed.matrix}} from a \code{.bed} file, and either
  a \code{.rds} file or a \code{.bim} and a \code{.fam} file.
}
\usage{
read.bed.matrix(basename, bed = paste(basename, ".bed", sep=""), 
                          fam = paste(basename, ".fam", sep=""),
                          bim = paste(basename, ".bim", sep=""), 
                          rds = paste(basename, ".rds", sep=""),
                          verbose = getOption("gaston.verbose",TRUE))
}

\arguments{
  \item{basename}{ Basename of all files }
  \item{bed}{ Name of the \code{.bed} file }
  \item{fam}{ Name of the \code{.fam} file }
  \item{bim}{ Name of the \code{.bim} file }
  \item{rds}{ Name of the \code{.rds} file (ignored if \code{NULL}) }
  \item{verbose}{ If \code{TRUE}, display information on the function actions }
}

\details{
If the \code{.rds} file exists, the \code{.fam} and \code{.bim} files will be 
ignored. To ignore an existing \code{.rds} file, set \code{rds = NULL}.

The \code{.bed}, \code{.fam} and \code{.bim} files follow the PLINK specifications
(\url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed}).

If the \code{.bed} file does not exist, and \code{basename} ends by \code{".bed"},
the function will try to generate a new basename by trimming the extension out. This
allows to write \code{read.bed.matrix("file.bed")} instead of \code{read.bed.matrix("file")}.
}
\value{ A \code{\link{bed.matrix}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{write.bed.matrix}} }

\examples{
# Read RDS and bed files
x <- read.bed.matrix( system.file("extdata", "LCT.bed", package="gaston") )
x
}
