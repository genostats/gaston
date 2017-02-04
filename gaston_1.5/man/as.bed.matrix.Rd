\name{as.bed.matrix}
\alias{as.bed.matrix}

\title{ Creation of a bed.matrix }

\description{ Creates a bed.matrix using a numeric matrix and two data frame for ped / snps slots }
\usage{ 
as.bed.matrix(x, fam, bim)
}

\arguments{
  \item{x}{ A numeric matrix }
  \item{fam}{ (Optionnal) A data frame (the contents of a .fam file) }
  \item{bim}{ (Optionnal) A data frame (the contents of a .bim file) }
}

\details{
The data frame \code{fam} should have columns named "famid", "id", "father", "mother", "sex" and "pheno".
The data frame \code{bim} should have columns named "chr", "id", "dist", "pos", "A1" and "A2". 
}

\value{
A bed.matrix condensing all three arguments.
}

\seealso{ \code{\link{bed.matrix-class}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\examples{
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
x
}

\keyword{ Association Test }
