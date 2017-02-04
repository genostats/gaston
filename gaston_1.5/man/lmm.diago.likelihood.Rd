\name{lmm.diago.likelihood}
\alias{lmm.diago.likelihood}

\title{ Likelihood of a linear mixed model }

\description{
  Compute the Restricted Likelihood of a linear mixed model, using the "diagonalization trick".
}
\usage{ lmm.diago.likelihood(tau, s2, h2, Y, X, eigenK, p = 0) }

\arguments{
  \item{tau}{ Value(s) of model parameter (see Details) }
  \item{s2}{  Value(s) of model parameter (see Details) }
  \item{h2}{  Value(s) of heritability (see Details) }
  \item{Y}{ Phenotype vector }
  \item{X}{ Covariable matrix }
  \item{eigenK}{ Eigen decomposition of \eqn{K} (a positive symmetric matrix)}
  \item{p}{ Number of Principal Components included in the mixed model with fixed effect }
}

\details{

  Compute the Restricted Likelihood under the linear mixed model
  \deqn{ Y = (X|PC)\beta + \omega + \varepsilon }{ Y = (X|PC) beta + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and 
  \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.

  The matrix \eqn{K} is given through its eigen decomposition, as produced by \code{eigenK = eigen(K, symmetric = TRUE)}.
  The matrix \eqn{(X|PC)} is the concatenation of the covariable matrix \eqn{X} and 
  of the first \eqn{p} eigenvectors of \eqn{K}, included in the model with fixed effects.
  
  If both \code{tau} and \code{s2} (for \eqn{\sigma^2}{sigma^2}) are provided, the function computes the
  likelihood for these values of the parameters; if these parameters are vectors of length \eqn{> 1},
  then a matrix of likelihood values is computed.

  If \code{h2} is provided, the function computes \eqn{\tau}{tau} and \eqn{\sigma^2}{s2} which
  maximizes the likelihood under the constraint \eqn{ {\tau \over \tau + \sigma^2 } = h^2 }{tau/(tau + s2) = h2},
  and outputs these values as well as the likelihood value at this point.
}

\value{
  If \code{tau} and \code{s2} are provided, the corresponding likelihood values.

  If \code{tau} or \code{s2} are missing, and \code{h2} is provided, a named list with members
    \item{ tau }{ Corresponding values of \eqn{\tau}{tau} }
    \item{ sigma2 }{ Corresponding values of \eqn{\sigma^2}{sigma^2}}
    \item{ likelihood }{ Corresponding likelihood values }
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{lmm.diago}}, \code{\link{lmm.aireml}} }

\examples{
# Load data
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
TAU <- seq(0.5,1.5,length=30)
S2 <- seq(1,3,length=30)
lik1 <- lmm.diago.likelihood(tau = TAU, s2 = S2, Y = y, eigenK = eiK)

H2 <- seq(0,1,length=51)
lik2 <- lmm.diago.likelihood(h2 = H2, Y = y, eigenK = eiK)

# Plotting
par(mfrow=c(1,2))
lik.contour(TAU, S2, lik1, heat = TRUE, xlab = "tau", ylab = "sigma^2")
lines(lik2$tau, lik2$sigma2)
plot(H2, exp(lik2$likelihood), type="l", xlab="h^2", ylab = "likelihood")
}

\keyword{ Eigen decomposition }
\keyword{ Likelihood }
