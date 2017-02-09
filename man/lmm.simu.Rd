\name{lmm.simu}
\alias{lmm.simu}

\title{ Linear mixed model data simulation}

\description{
  Simulate data under a linear mixed model, using the eigen decomposition
  of the variance matrix.
}
\usage{ lmm.simu(tau, sigma2, K, eigenK = eigen(K), X, beta) }

\arguments{
  \item{tau}{ Model parameter }
  \item{sigma2}{ Model parameter }
  \item{K}{ (Optional) A positive symmetric matrix \eqn{K} }
  \item{eigenK}{ Eigen decomposition of \eqn{K} }
  \item{X}{ Covariable matrix }
  \item{beta}{ Fixed effect vector of covariables }
}

\details{
  The data are simulated under the following linear mixed model :
  \deqn{ Y = X\beta + \omega + \varepsilon }{ Y = X beta + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and
  \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.

  The simulation uses \eqn{K} only through its eigen decomposition; the parameter
  \code{K} is therefore optional.
}

\value{
  A named list with two members:
    \item{y}{ Simulated value of \eqn{Y} }
    \item{omega}{Simulated value of \eqn{\omega}{omega}  }
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{random.pm}} }

\examples{
# generate a random positive matrix 
set.seed(1)
R <- random.pm(503)

# simulate data with a "polygenic component" 
y <-  lmm.simu(0.3, 1, eigenK = R$eigen)
str(y)
}

\keyword{ Simulations }
