\name{score.fixed.linear/score.fixed.logistic}
\alias{score.fixed.linear}
\alias{score.fixed.logistic}

\title{ Score Test for Covariates with Fixed Effects in Linear or Logistic Mixed Model }

\description{ Score Test for association between covariates and phenotype.}
\usage{ 
score.fixed.linear(x, Y, X = matrix(1, length(Y)), K, ...)
score.fixed.logistic(x, Y, X = matrix(1, length(Y)), K, ...)
}

\arguments{
  \item{x}{ A matrix of covariates }
  \item{Y}{ The phenotype vector }
  \item{X}{ A covariable matrix. The default is a column vector of ones, to include an intercept in the model }
  \item{K}{ A positive definite matrix or a \code{list} of such matrices }
  \item{...}{ Optional arguments used to fit null model in \code{lmm.aireml} or \code{logistic.mm.aireml} function. }
}

\details{
  The function \code{score.fixed.linear} considers the linear mixed model
  \deqn{ Y = X\alpha + x\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y = X alpha + x beta + omega_1 + ... + omega_k + varepsilon }
  whereas the \code{score.fixed.logistic} function considers the following logistic model
  \deqn{ \mbox{logit}(P[Y=1|X,x,\omega_1,\ldots,\omega_k])  = X\alpha + x\beta + \omega_1 + \ldots + \omega_k}{logit(P[Y=1|X,x,omega_1,...,omega_k])  = Xalpha + x beta + omega_1 + ... + omega_k}
  with \eqn{ \omega_j \sim N(0,\tau_j K_j) }{omega_j ~ N(0, tau_j K_j)} where \eqn{K_j} are Genetic Relationship Matrix (GRM), \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}
  and fixed effects \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
 
  The two functions give score test for
  \eqn{H_0}{H0} : \eqn{\beta=0}{beta=0} vs \eqn{H_1}{H1} : \eqn{\beta\neq 0}{beta !=0}.  
  In this aim, all parameters under null model are estimated with \code{lmm.aireml} or \code{logistic.mm.aireml}.
}

\value{
  A named list of values:
  \item{score}{ Estimated score }
  \item{p}{ The corresponding p-value }
  \item{log.p}{ The logarithm of corresponding p-value }
}

\seealso{ \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p"

# Calculate GRM et its eigen decomposition
k <- GRM(x)
eig <- eigen(k)
eig$values <- round(eig$values, 5)

# generate covariate matrix
set.seed(1)
X <- cbind( rbinom(nrow(x), 1, prob=1/2), rnorm(nrow(x)) )


# simulate quantitative phenotype with polygenic component and covariate effects
y <- X \%*\% c(-1,0.5) + lmm.simu(0.3,1,eigenK=eig)$y

t <- score.fixed.linear(X, y, K=k, verbose=FALSE)
str(t)


# simulate binary phenotype with polygenic component and covariate effects
mu <- X \%*\% c(-1,0.5) + lmm.simu(1, 0, eigenK=eig)$y
pi <- 1/(1+exp(-mu))
y <- 1*( runif(length(pi))<pi )

tt <- score.fixed.logistic(X, y, K=k, verbose=FALSE)
str(tt)
}

\keyword{ Score Test }
