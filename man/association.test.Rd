\name{association.test}
\alias{association.test}

\title{ Association Test }

\description{ Association tests between phenotype and SNPs. }
\usage{ 
association.test(x, Y = x@ped$pheno, X = matrix(1, nrow(x)),
                 method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                 test = c("score", "wald", "lrt"), K, eigenK, beg = 1, 
                 end = ncol(x), p = 0, tol = .Machine$double.eps^0.25, ...)
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{Y}{ The phenotype vector. Default is the column (\code{pheno}) of \code{x@ped} }
  \item{X}{ A covariable matrix. The default is a column vector of ones, to include an intercept in the model }
  \item{method}{ Method to use: \code{"lm"} for (generalized) linear model, and \code{"lmm"} for (generalized) linear mixed model }
  \item{response}{ Is \code{"Y"} a quantitative or a binary phenotype?}
  \item{test}{ Which test to use. For binary phenotypes, \code{test = "score"} is mandatory}
  \item{K}{ A Genetic Relationship Matrix (as produced by \code{\link{GRM}}), or a list of such matrices. For \code{test = "score"}. }
  \item{eigenK}{ Eigen decomposition of the Genetic Relationship Matrix (as produced by the function \code{eigen}). 
                 For \code{test = "wald"} or \code{"lrt"}. }
  \item{beg}{ Index of the first SNP tested for association }
  \item{end}{ Index of the last SNP tested for association }
  \item{p}{ Number of Principal Components to include in the model with fixed effect (for \code{test = "wald"} or \code{"lrt"}) }
  \item{tol}{ Parameter for the likelihood maximization (as in \code{optimize}) }
  \item{...}{ Additional parameters for \code{\link{lmm.aireml}} or \code{\link{logistic.mm.aireml}} (if \code{test = "score"}). }
}

\details{
  Tests the association between the phenotype and requested SNPs in \code{x}.

  If \code{method = "lm"} and \code{response = "quantitative"} are used, a simple linear regression 
  is performed to test each SNP in the considered interval. Precisely, the following model is
  considered for each SNP,
  \deqn{ Y = (X|PC)\alpha + G\beta + \varepsilon }{ Y = (X|PC) alpha + G beta + epsilon }
  with \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)},
  \eqn{G} the genotype vector of the SNP, 
  \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.
  A Wald test is performed, independently of the value of \code{test}.

  If\code{method = "lm"} and \code{response = "binary"}, a similar model is used for a logistic 
  regression (Wald test).
 
  If \code{method = "lmm"} and \code{response = "quantitative"}, the following model in considered for each SNP
  \deqn{ Y = (X|PC)\alpha + G\beta + \omega + \varepsilon }{ Y = (X|PC) alpha + G beta + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
  \eqn{G} is the genotype vector of the SNP, \eqn{K} is a Genetic Relationship Matrix (GRM)
  \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.

  If \code{test = "score"}, all parameters are estimated with the same procedure as in 
  \code{\link{lmm.aireml}} and the argument \code{K} is used to specify the GRM matrix (or a list of GRM
  matrices for inclusion of several random effects in the model). If \code{p} is positive, the paramater \code{eigenK}
  needs to be given as well.
  For Wald and LRT tests the procedure used is the same as in \code{\link{lmm.diago}} and \code{eigenK} is used to 
  specify the GRM matrix.

  If \code{method = "lmm"} and \code{response = "binary"}, the following model in considered for each SNP
  \deqn{ \mbox{logit}(P[Y=1| X, G, \omega])  = X\alpha + G\beta + \omega}{logit P(Y=1|X,G,omega)  = X alpha + G beta + omega}
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)}. 
  \eqn{G} is the genotype vector of the SNP, \eqn{K}{K} is a Genetic Relationship Matrix (GRM),
  \eqn{X} the covariable matrix. A score test is performed, independently of the value of \code{test}.
  All parameters under null model are estimated with the same procedure as in \code{\link{logistic.mm.aireml}}.
  In case of convergence problems of the null problem, the user can try several starting values (in particular
  with parameter \code{tau}, trying e.g. \code{tau = 0.1} or another value).
  It is possible to give a list of matrices in parameter \code{K} for inclusion of several random effects in the model.
  If \code{p} is positive, the paramater \code{eigenK} needs to be given as well.

  Note: this function is not multithreaded. Wald test with Linear Mixed Models are computationally intensive, 
  to run a GWAS with such tests consider using \code{association.test.parallel} in package \code{gaston.utils} 
  (on github). Association tests with dosages can be done with \code{association.test.dosage} and 
  \code{association.test.dosage.parallel} in the same package.
}

\value{
  A data frame, giving for each considered SNP, its position, id, alleles, and 
  some of the following columns depending on the values of \code{method} and \code{test}:
  \item{score}{Score statistic for each SNP}
  \item{h2}{Estimated value of \eqn{\tau \over {\tau + \sigma^2}}{tau/(tau + sigma^2)}}
  \item{beta}{Estimated value of \eqn{\beta}{beta}}
  \item{sd}{Estimated standard deviation of the \eqn{\beta}{beta} estimation}
  \item{LRT}{Value of the Likelihood Ratio Test}
  \item{p}{The corresponding p-value}
}

\seealso{ \code{\link{qqplot.pvalues}}, \code{\link{manhattan}},  \code{\link{lmm.diago}},  
          \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}}, \code{\link[stats:optimize]{optimize}} }

\examples{
\donttest{
# Load data
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
standardize(x) <- "p"

# simulate quantitative phenotype with effect of SNP #631
set.seed(1)
y <- x \%*\% c(rep(0,630),0.5,rep(0,ncol(x)-631)) + rnorm(nrow(x))

# association test with linear model 
test <- association.test(x, y, method="lm", response = "quanti")

# a p-values qq plot
qqplot.pvalues(test)

# a small Manhattan plot 
# hihlighting the link between p-values and LD with SNP #631
lds <- LD(x, 631, c(1,ncol(x)))
manhattan(test, col = rgb(lds,0,0), pch = 20)

# use y to simulate a binary phenotype
y1 <- as.numeric(y > mean(y))

# logistic regression
t_binary <- association.test(x, y1, method = "lm", response = "binary")
# another small Manhattan plot
manhattan(t_binary, col = rgb(lds,0,0), pch = 20)

}}

\keyword{ Association Test }
