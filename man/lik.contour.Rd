\name{lik.contour}
\alias{lik.contour}

\title{ Contour plot for two parameters likelihood }

\description{ Create a contour plot (superimposed with a heat map) }
\usage{ lik.contour(x, y, z, levels = NULL, nlevels = 11, heat = TRUE, col.heat = NULL, ...) }

\arguments{
  \item{x, y, z}{ As in \code{contour} }
  \item{levels}{ As in \code{contour}. If \code{NULL}, the function computes appropriate levels. }
  \item{nlevels}{ As in \code{contour} }
  \item{heat}{ If \code{TRUE}, a heat map is superimposed to the contour plot }
  \item{col.heat}{ Vector of heat colors}
  \item{\dots}{ Additional arguments to \code{image} and \code{contour}}
}
\details{ This function is a wrapper for \code{contour}, with a different method to compute
a default value for levels. If \code{heat = TRUE}, a heatmap produced by \code{image} is added to the plot.
See \code{\link{contour}} for details on parameters. }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{lmm.diago.likelihood}}, \code{\link[graphics:contour]{contour}}, \code{\link[graphics:image]{image}} }

\examples{
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)

# Compute Genetic Relationship Matrix
K <- GRM(x)

# eigen decomposition of K
eiK <- eigen(K)

# simulate a phenotype
set.seed(1)
y <- 1 + lmm.simu(tau = 1, sigma2 = 2, eigenK = eiK)$y

# Likelihood
TAU <- seq(0.5,2.5,length=30)
S2 <- seq(1,3,length=30)
lik1 <- lmm.diago.likelihood(tau = TAU, s2 = S2, Y = y, eigenK = eiK)
lik.contour(TAU, S2, lik1, heat = TRUE, xlab = "tau", ylab = "sigma^2")
}

\keyword{ Heat map }

