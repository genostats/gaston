\name{score.variance.linear/score.variance.logistic}
\alias{score.variance.linear}
\alias{score.variance.logistic}

\title{ Variance Component Test in Linear or Logistic Mixed Model}

\description{ Test if a variance component is significaly different from 0 using score test in a Linear or Logistic Mixed Model. }
\usage{ 
score.variance.linear(K0, Y, X = matrix(1, length(Y)), K, acc_davies=1e-10, ...)
score.variance.logistic(K0, Y, X = matrix(1, length(Y)), K, acc_davies=1e-10, ...)
}

\arguments{
  \item{K0}{ A positive definite matrix }
  \item{Y}{ The phenotype vector }
  \item{X}{ A covariate matrix. The default is a column vector of ones, to include an intercept in the model }
  \item{K}{ A positive definite matrix or a \code{list} of such matrices }
  \item{acc_davies}{ Accuracy in Davies method used to compute p-value }
  \item{...}{ Optional arguments used to fit null model with \code{lmm.aireml} of \code{logistic.mm.aireml} function. }
}

\details{
  In \code{score.variance.linear}, we consider the linear mixed model
  \deqn{ Y = X\alpha + \gamma + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y = X alpha + gamma + omega_1 + ... + omega_k + varepsilon }
  or, in \code{score.variance.logistic}, we consider the following logistic model
  \deqn{ \mbox{logit}(P[Y=1|X,x,\omega_1,\ldots,\omega_k]) = X\alpha + \gamma + \omega_1 + \ldots + \omega_k}{logit(P[Y=1|X,x,omega_1,...,omega_k]) = X alpha + gamma + omega_1 + ... + omega_k }
  with \eqn{ \gamma\sim N(0,\kappa K_0)\gamma}{gamma~N(0, kappa K_0)}, \eqn{ \omega_j \sim N(0,\tau_j K_j) }{omega_j ~ N(0, tau_j K_j)},
  \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
  \eqn{K_0} and \eqn{K_j} are Genetic Relationship Matrix (GRM).  
  
  \code{score.variance.linear } and \code{score.variance.logistic} functions permit to test
  \deqn{H_0 : \kappa=0 \mbox{ vs } H_1 : \kappa>0}{H_0 : kappa=0 vs H_1 : kappa>0}
  with, for linear mixed model, the score
  \deqn{ Q = Y'P_OK_0P_0Y/2 }{ Q = Y'P_OK_0P_0Y/2 }
  or, for logistic mixed model, the score
  \deqn{ Q = (Y-\pi_0)'K_0(Y-\pi_0)/2 }{ Q = (Y-pi_0)'K_0(Y-pi_0)/2 }
  where \eqn{P_0} is the last matrix \eqn{P} computed in the optimization process for null model and \eqn{\pi_0}{pi_0} the vector of fitted values under null logistic model.
  
  The associated p-value is computed with Davies method.
  
  
  In this aim, all parameters under null model are estimated with \code{lmm.aireml} or \code{logistic.mm.aireml}.
  The p-value corresponding to the estimated score is computed using Davies method implemented in 'CompQuadForm' R package.
}

\value{
  A named list of values:
  \item{score}{ Estimated score }
  \item{p}{ The corresponding p-value }
}

\seealso{ \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}} }

\references{Davies R.B. (1980) \emph{Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables},
                Journal of the Royal Statistical Society. Series C (Applied Statistics), \bold{323-333}}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p"

# Calculate GRM et its eigen decomposition
K0 <- GRM(x)
eig <- eigen(K0)
eig$values <- round(eig$values, 5)

# generate an other positive matrix (to play the role of the second GRM)
set.seed(1)
R <- random.pm(nrow(x))


# simulate quantitative phenotype with two polygenic components
y <- lmm.simu(0.1,1,eigenK=eig)$y + lmm.simu(0.2,0,eigenK=R$eigen)$y

t <- score.variance.linear(K0, y, K=R$K, verbose=FALSE)
str(t)


# simulate binary phenotype with two polygenic components
mu <- lmm.simu(0.1,0.5,eigenK=eig)$y + lmm.simu(0.2,0,eigenK=R$eigen)$y
pi <- 1/(1+exp(-mu))
y <- 1*(runif(length(pi))<pi)

tt <- score.variance.logistic(K0, y, K=R$K, verbose=FALSE)
str(tt)
}

\keyword{ Variance Component Test }
