\name{logistic.mm.aireml}
\alias{logistic.mm.aireml}

\title{ Logistic mixed model fitting with Penalized Quasi-Likelihood / AIREML }

\description{
Estimate the parameters of a logistic linear mixed model using the Penalized
Quasi-Likelihood with an AIREML step for the linear model.
}
\usage{
logistic.mm.aireml(Y, X = matrix(1, nrow = length(Y)), K, 
                   min_tau, tau, beta, constraint = TRUE, max.iter = 50L, eps = 1e-5,
                   verbose = getOption("gaston.verbose",TRUE), get.P = FALSE, EM = FALSE)
}
\arguments{
  \item{Y}{ Binary phenotype vector }
  \item{X}{ Covariable matrix. By default, a column of ones to include an intercept in the model }
  \item{K}{ A positive definite matrix or a \code{list} of such matrices }
  \item{min_tau}{ Minimal value for model parameter \eqn{\tau}{tau} (if missing, will be set to \eqn{10^{-6}}{1e-6})}
  \item{tau}{ (Optional) Optimization starting point for variance component(s) \code{tau} }
  \item{beta}{ (Optional) Optimization starting point for fixed effect(s) \code{beta} }
  \item{constraint}{ If \code{TRUE}, the model parameters respect the contraints given by \code{min_tau} }
  \item{max.iter}{ Maximum number of iterations }
  \item{eps}{ The algorithm stops when the gradient norm is lower than this parameter }
  \item{verbose}{ If \code{TRUE}, display information on the algorithm progress }
  \item{get.P}{If \code{TRUE}, the function sends back the last matrix \eqn{P} computed in the optimization process}
  \item{EM}{If \code{TRUE}, the AIREML step is replaced by an EM step}
}

\details{
  Estimate the parameters of the following logistic mixed model:
  \deqn{ \mbox{logit}(P[Y=1|X,\omega_1,\ldots,\omega_k])  = X\beta + \omega_1 + \ldots + \omega_k}
  {logit P(Y=1|X,omega_1,...,omega_k)  = X beta  + omega_1 + ... + omega_k}
  with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k }.
  
  The estimation is based on the Penalized Quasi-Likelihood with an AIREML step for the linear model
  (the algorithm is similar to the algorithm described in Chen et al 2016). If \code{EM = TRUE}
  the AIREML step is replaced by an EM step. In this case the convergence will be much slower, 
  you're advised to use a large value of \code{max.iter}.

  The variance matrices \eqn{K_1}{K_1}, ..., \eqn{K_k}{K_k}, are specified through the parameter \code{K}.

  After convergence, the function also compute Best Linear Unbiased Predictors (BLUPs) for 
  \eqn{\beta}{beta} and \eqn{\omega}{omega}.
}

\value{
  A named list with members:
    \item{tau}{ Estimate(s) of the model parameter(s) \eqn{\tau_1, \dots, \tau_k}{tau_1, ..., tau_k} }
    \item{niter}{ Number of iterations done }
    \item{P}{ Last computed value of matrix P (see reference) }
    \item{BLUP_omega}{ BLUPs of random effects }
    \item{BLUP_beta}{ BLUPs of fixed effects \eqn{\beta}{beta}}
    \item{varbeta}{ Variance matrix for \eqn{\beta}{beta} estimates }
  If \code{get.P = TRUE}, there is an additional member:
    \item{P}{The last matrix \eqn{P} computed in the AIREML step}
}

\references{Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995), \emph{Average information REML: an efficient algorithm for variance
parameter estimation in linear mixed models}, Biometrics, \bold{1440-1450}

Chen, Han et al. (2016), \emph{Control for Population Structure and Relatedness for Binary Traits in Genetic Association Studies via
Logistic Mixed Models}, The American Journal of Human Genetics, \bold{653--666}}

\seealso{ \code{\link{lmm.aireml}}, \code{\link{lmm.diago}}, \code{\link{lmm.simu}}}

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)

# Compute Genetic Relationship Matrix
standardize(x) <- "p"
K <- GRM(x)

# Simulate a quantitative genotype under the LMM
set.seed(1)
mu <- 1 + x \%*\% rnorm(ncol(x), sd = 2)/sqrt(ncol(x))
pi <- 1/(1+exp(-mu))
y <- 1*( runif(length(pi))<pi )

# Estimates
estimates <- logistic.mm.aireml(y, K = K, verbose = FALSE)
str(estimates)
}

\keyword{ Penalized Quasi-Likelihood }
\keyword{ Average Information Restricted Maximum Likelihood (AIREML) }
\keyword{ Genetic effect }
