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
The \code{.bed}, \code{.fam} and \code{.bim} files follow the PLINK specifications
(\url{http://zzz.bwh.harvard.edu/plink/binary.shtml}).

If a \code{.rds} file exists (created by \code{write.bed.matrix}), 
the \code{.fam} and \code{.bim} files will be ignored. 
To ignore an existing \code{.rds} file, set \code{rds = NULL}.

If the \code{.bed} file does not exist, and \code{basename} ends by \code{".bed"},
the function will try to generate a new basename by trimming the extension out. This
allows to write \code{read.bed.matrix("file.bed")} instead of \code{read.bed.matrix("file")}.

If the option \code{gaston.auto.set.stats} is set to \code{TRUE} (the default), 
the function \code{set.stats} will be called before returning the \code{bed.matrix},
unless a \code{.rds} file is present: in this case, the \code{bed.matrix} obtained
is identical to the \code{bed.matrix} saved with \code{write.bed.matrix}.
}
\value{ A \code{\link{bed.matrix}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{write.bed.matrix}}, \code{\link{set.stats}} }

\examples{
# Read RDS and bed files
x <- read.bed.matrix( system.file("extdata", "LCT.bed", package="gaston") )
x
}
