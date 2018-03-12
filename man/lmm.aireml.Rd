\name{lmm.aireml}
\alias{lmm.aireml}

\title{ Linear mixed model fitting with AIREML }

\description{
Estimate the parameters of a linear mixed model, 
using Average Information Restricted Maximum Likelihood (AIREML) algorithm.
}
\usage{
lmm.aireml(Y, X = matrix(1, nrow = length(Y)), K,
           EMsteps = 0L, EMsteps_fail = 1L, EM_alpha = 1,
           min_tau, min_s2 = 1e-06, theta, constraint = TRUE, max_iter = 50L,
           eps = 1e-05, verbose = getOption("gaston.verbose", TRUE),
           contrast = FALSE, get.P = FALSE) }

\arguments{
  \item{Y}{ Phenotype vector }
  \item{X}{ Covariable matrix. By default, a column of ones to include an intercept in the model }
  \item{K}{ A positive definite matrix or a \code{list} of such matrices }
  \item{EMsteps}{ Number of EM steps ran prior the AIREML }
  \item{EMsteps_fail}{ Number of EM steps performed when the AIREML algorithm fail to improve the likelihood value }
  \item{EM_alpha}{ Tweaking parameter for the EM (see Details) }
  \item{min_tau}{ Minimal value for model parameter \eqn{\tau}{tau} (if missing, will be set to \eqn{10^{-6}}{1e-6})}
  \item{min_s2}{ Minimal value for model parameter \eqn{\sigma^2}{sigma^2} }
  \item{theta}{ (Optional) Optimization starting point \code{theta = c(sigma^2, tau)} }
  \item{constraint}{ If \code{TRUE}, the model parameters respect the contraints given by \code{min_tau} and \code{min_s2} }
  \item{max_iter}{ Maximum number of iterations }
  \item{eps}{ The algorithm stops when the gradient norm is lower than this parameter }
  \item{verbose}{ If \code{TRUE}, display information on the algorithm progress }
  \item{contrast}{ If \code{TRUE}, use a contrast matrix to compute the Restricted Likelihood (usually slower)}
  \item{get.P}{If \code{TRUE}, the function sends back the last matrix \eqn{P} computed in the optimization process}
}

\details{
  Estimate the parameters of the following linear mixed model, using AIREML algorithm:
  \deqn{ Y = X\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y =X beta + omega_1 + ... + omega_k + epsilon }
  with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k } and 
  \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.

  The variance matrices \eqn{K_1}{K_1}, ..., \eqn{K_k}{K_k}, are specified through the parameter \code{K}.

  If \code{EMsteps} is positive, the function will use this number of EM steps to compute a better starting point
  for the AIREML algorithm. Setting \code{EMsteps} to a value higher than \code{max_iter} leads to an EM optimization.
  It can happen that after an AIREML step, the likelihood did not increase: if this 
  happens, the functions falls back to \code{EMsteps_fail} EM steps. The parameter \code{EM_alpha} can be set to 
  a value higher than \eqn{1} to attempt to accelerate EM convergence; this could also result in uncontrolled
  behaviour and should be used with care.

  After convergence, the function also compute Best Linear Unbiased Predictors (BLUPs) for 
  \eqn{\beta}{beta} and \eqn{\omega}{omega}, and an 
  estimation of the participation of the fixed effects to the variance of \eqn{Y}.
}

\value{
  A named list with members:
    \item{sigma2}{ Estimate of the model parameter \eqn{\sigma^2}{sigma^2} }
    \item{tau}{ Estimate(s) of the model parameter(s) \eqn{\tau_1, \dots, \tau_k}{tau_1, ..., tau_k} }
    \item{logL}{ Value of log-likelihood }
    \item{logL0}{ Value of log-likelihood under the null model (without random effect) }
    \item{niter}{ Number of iterations done }
    \item{norm_grad}{ Last computed gradient's norm }
    \item{P}{ Last computed value of matrix P (see reference) }
    \item{Py}{ Last computed value of vector Py (see reference) }
    \item{BLUP_omega}{ BLUPs of random effects }
    \item{BLUP_beta}{ BLUPs of fixed effects \eqn{\beta}{beta}}
    \item{varbeta}{ Variance matrix for \eqn{\beta}{beta} estimates }
    \item{varXbeta}{ Participation of fixed effects to variance of Y }
  If \code{get.P = TRUE}, there is an additional member:
    \item{P}{The last matrix \eqn{P} computed in the AIREML step}
}

\references{Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995), \emph{Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models}, Biometrics, \bold{1440-1450}}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{\code{\link{lmm.diago}}, \code{\link{logistic.mm.aireml}}, \code{\link{lmm.simu}}}

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)

# Compute Genetic Relationship Matrix
standardize(x) <- "p"
K <- GRM(x)

# Simulate a quantitative genotype under the LMM
set.seed(1)
y <- 1 + x \%*\% rnorm(ncol(x), sd = 1)/sqrt(ncol(x)) + rnorm(nrow(x), sd = sqrt(2))

# Estimates
estimates <- lmm.aireml(y, K = K, verbose = FALSE)
str(estimates)
}

\keyword{ Average Information Restricted Maximum Likelihood (AIREML) }
\keyword{ Genetic effect }
