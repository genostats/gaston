\name{association.test}
\alias{association.test}

\title{ Association Test }

\description{ Association tests between phenotype and SNPs. }
\usage{ 
association.test(x, Y = x@ped$pheno, X = matrix(1, nrow(x)),
                 method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                 test = c("score", "wald", "lrt"),
                 K, eigenK, beg = 1, end = ncol(x), p = 0,
                 tol = .Machine$double.eps^0.25, multithreaded = FALSE, ...)
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
  \item{multithreaded}{ \code{Logical}. If \code{TRUE}, several threads are used. }
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

  There is no test available yet for \code{method = "lm"} and \code{response = "binary"}. You can consider
  treating the binary phenotype as quantitative.
 
  If \code{method = "lmm"} and \code{response = "quantitative"}, the following model in considered for each SNP
  \deqn{ Y = (X|PC)\alpha + G\beta + \omega + \varepsilon }{ Y = (X|PC) alpha + G beta + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
  \eqn{G} is the genotype vector of the SNP, \eqn{K} is a Genetic Relationship Matrix (GRM)
  \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.

  If \code{test = "score"}, all parameters are estimated with the same procedure as in 
  \code{\link{lmm.aireml}} and the argument \code{K} is used to specify the GRM matrix (or a list of GRM
  matrices for inclusion of several random effects in the model). The argument \code{p} is ignored.
  For Wald and LRT tests the procedure used is the same as in \code{\link{lmm.diago}} and \code{eigenK} is used to 
  specify the GRM matrix.

  If \code{method = "lmm"} and \code{response = "binary"}, the following model in considered for each SNP
  \deqn{ \mbox{logit}(P[Y=1| X, G, \omega])  = X\alpha + G\beta + \omega}{logit P(Y=1|X,G,omega)  = X alpha + G beta + omega}
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)}. 
  \eqn{G} is the genotype vector of the SNP, \eqn{K}{K} is a Genetic Relationship Matrix (GRM),
  \eqn{X} the covariable matrix. A score test is performed, independantly of the value of \code{test}.
  All parameters under null model are estimated with the same procedure as in \code{\link{logistic.mm.aireml}}.
  It is possible to give a list in argument \code{K} for inclusion of several random effects in the model.
  The argument \code{p} is ignored.

  Note: in the present implementation, using several threads is possible for Wald and LRT tests, but is
  likely to be slower than using only one, hence the \code{multithreaded} parameter. This will change in the future.
}

\value{
  A data frame, giving for each considered SNP, some of the following columns depending on the 
  values of \code{method} and \code{test}:
  \item{score}{Score statistic for each SNP}
  \item{h2}{Estimated value of \eqn{\tau \over {\tau + \sigma^2}}{tau/(tau + sigma^2)}}
  \item{beta}{Estimated value of \eqn{\beta}{beta}}
  \item{sd}{Estimated standard deviation of the \eqn{\beta}{beta} estimation}
  \item{LRT}{Value of the Likelihood Ratio Test}
  \item{p}{The corresponding p-value}
}

\seealso{ \code{\link{lmm.diago}},  \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}}, \code{\link[stats:optimize]{optimize}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\examples{
\donttest{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
standardize(x) <- "p"

# generate a random positive matrix (to play the role of the GRM)
set.seed(1)
R <- random.pm(nrow(x))

# simulate quantitative phenotype with effect of SNP #351 and a polygenic component 
y <- x \%*\% c(rep(0,350),0.25,rep(0,ncol(x)-351)) + lmm.simu(0.3,1,eigenK=R$eigen)$y

# association test with mixed model (default is quantitative response, score test)
t_score <- association.test(x, y, K = R$K, method = "lmm")
t_wald <- association.test(x, y, eigenK = R$eigen, method = "lmm", test = "wald")

# comparaison des p-valeurs obtenues par les deux tests
plot( t_score$p, t_wald$p, log = "xy", xlab = "score", ylab = "wald")

# (mini) Manhattan plot
plot(-log10(t_score$p), xlab="SNP index", ylab = "-log(p)")

# link between p-values and LD with SNP #351
lds <- LD(x, 351, c(1,ncol(x)))
plot(lds, -log10(t_score$p), xlab="r^2", ylab="-log(p)")

# use y to simulate a binary phenotype
y1 <- as.numeric(y > mean(y))
t_binary <- association.test(x, y1, K = R$K, response = "binary")
# (mini) Manhattan plot
plot(-log10(t_binary$p), xlab="SNP index", ylab = "-log(p)")

}}

\keyword{ Association Test }
