#' Association Test
#' 
#' @description  Association tests between phenotype and SNPs.
#' 
#' @param x  A \code{\link{bed.matrix}}
#' @param Y  The phenotype vector. Default is the column (\code{pheno}) of \code{x@ped}
#' @param X  A covariable matrix. The default is a column vector of ones, to include an intercept in the model
#' @param method  Method to use: \code{"lm"} for (generalized) linear model, and \code{"lmm"} for (generalized) linear mixed model
#' @param response  Is \code{"Y"} a quantitative or a binary phenotype?
#' @param test  Which test to use. For binary phenotypes, \code{test = "score"} is mandatory
#' @param K  A Genetic Relationship Matrix (as produced by \code{\link{GRM}}), or a list of such matrices. For \code{test = "score"}.
#' @param eigenK  Eigen decomposition of the Genetic Relationship Matrix (as produced by the function \code{eigen}).
#' For \code{test = "wald"} or \code{"lrt"}.
#' @param beg  Index of the first SNP tested for association
#' @param end  Index of the last SNP tested for association
#' @param p  Number of Principal Components to include in the model with fixed effect (for \code{test = "wald"} or \code{"lrt"})
#' @param tol  Parameter for the likelihood maximization (as in \code{optimize})
#' @param ...  Additional parameters for \code{\link{lmm.aireml}} or \code{\link{logistic.mm.aireml}} (if \code{test = "score"}).
#' 
#' 
#' @details
#' Tests the association between the phenotype and requested SNPs in \code{x}.
#' 
#' If \code{method = "lm"} and \code{response = "quantitative"} are used, a simple linear regression
#' is performed to test each SNP in the considered interval. Precisely, the following model is
#' considered for each SNP,
#' \deqn{ Y = (X|PC)\alpha + G\beta + \varepsilon }{ Y = (X|PC) alpha + G beta + epsilon }
#' with \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)},
#' \eqn{G} the genotype vector of the SNP,
#' \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.
#' A Wald test is performed, independently of the value of \code{test}.
#' 
#' If\code{method = "lm"} and \code{response = "binary"}, a similar model is used for a logistic
#' regression (Wald test).
#' 
#' If \code{method = "lmm"} and \code{response = "quantitative"}, the following model in considered for each SNP
#' \deqn{ Y = (X|PC)\alpha + G\beta + \omega + \varepsilon }{ Y = (X|PC) alpha + G beta + omega + epsilon }
#' with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
#' \eqn{G} is the genotype vector of the SNP, \eqn{K} is a Genetic Relationship Matrix (GRM)
#' \eqn{X} the covariates matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.
#' 
#' If \code{test = "score"}, all parameters are estimated with the same procedure as in
#' \code{\link{lmm.aireml}} and the argument \code{K} is used to specify the GRM matrix (or a list of GRM
#' matrices for inclusion of several random effects in the model). If \code{p} is positive, the paramater \code{eigenK}
#' needs to be given as well.
#' For Wald and LRT tests the procedure used is the same as in \code{\link{lmm.diago}} and \code{eigenK} is used to
#' specify the GRM matrix.
#' 
#' If \code{method = "lmm"} and \code{response = "binary"}, the following model in considered for each SNP
#' \deqn{ \mbox{logit}(P[Y=1| X, G, \omega])  = X\alpha + G\beta + \omega}{logit P(Y=1|X,G,omega)  = X alpha + G beta + omega}
#' with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)}.
#' \eqn{G} is the genotype vector of the SNP, \eqn{K}{K} is a Genetic Relationship Matrix (GRM),
#' \eqn{X} the covariable matrix. A score test is performed, independently of the value of \code{test}.
#' All parameters under null model are estimated with the same procedure as in \code{\link{logistic.mm.aireml}}.
#' In case of convergence problems of the null problem, the user can try several starting values (in particular
#' with parameter \code{tau}, trying e.g. \code{tau = 0.1} or another value).
#' It is possible to give a list of matrices in parameter \code{K} for inclusion of several random effects in the model.
#' If \code{p} is positive, the paramater \code{eigenK} needs to be given as well.
#' 
#' Note: this function is not multithreaded. Wald test with Linear Mixed Models are computationally intensive,
#' to run a GWAS with such tests consider using \code{association.test.parallel} in package \code{gaston.utils}
#' (on github). Association tests with dosages can be done with \code{association.test.dosage} and
#' \code{association.test.dosage.parallel} in the same package.
#' 
#' 
#' 
#' @return
#' A data frame, giving for each considered SNP, its position, id, alleles, and
#' some of the following columns depending on the values of \code{method} and \code{test}:
#' \item{score}{Score statistic for each SNP}
#' \item{h2}{Estimated value of \eqn{\tau \over {\tau + \sigma^2}}{tau/(tau + sigma^2)}}
#' \item{beta}{Estimated value of \eqn{\beta}{beta}}
#' \item{sd}{Estimated standard deviation of the \eqn{\beta}{beta} estimation}
#' \item{LRT}{Value of the Likelihood Ratio Test}
#' \item{p}{The corresponding p-value}
#' @seealso  \code{\link{qqplot.pvalues}}, \code{\link{manhattan}},  \code{\link{lmm.diago}},
#' \code{\link{lmm.aireml}}, \code{\link{logistic.mm.aireml}}, \code{\link[stats:optimize]{optimize}}
#' 
#' 
#' @keywords  Association Test
#' @examples
#' 
#' \donttest{
#' # Load data
#' data(TTN)
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' standardize(x) <- "p"
#' 
#' # simulate quantitative phenotype with effect of SNP #631
#' set.seed(1)
#' y <- x %*% c(rep(0,630),0.5,rep(0,ncol(x)-631)) + rnorm(nrow(x))
#' 
#' # association test with linear model 
#' test <- association.test(x, y, method="lm", response = "quanti")
#' 
#' # a p-values qq plot
#' qqplot.pvalues(test)
#' 
#' # a small Manhattan plot 
#' # hihlighting the link between p-values and LD with SNP #631
#' lds <- LD(x, 631, c(1,ncol(x)))
#' manhattan(test, col = rgb(lds,0,0), pch = 20)
#' 
#' # use y to simulate a binary phenotype
#' y1 <- as.numeric(y > mean(y))
#' 
#' # logistic regression
#' t_binary <- association.test(x, y1, method = "lm", response = "binary")
#' # another small Manhattan plot
#' manhattan(t_binary, col = rgb(lds,0,0), pch = 20)
#' 
#' }
#' @export association.test
association.test <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), 
                             method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                             test = c("score", "wald", "lrt"), 
                             K, eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, ...) {

  if(any(is.na(Y))) stop("This function does not accomodate missing phenotypes")
  
  response <- match.arg(response)
  if(response == "binary") {
    if(any(Y != 0 & Y != 1)) stop("Binary response should be 0 or 1")
    Y <- as.integer(Y)
  }

  if(beg < 1 | end > ncol(x)) stop("range too wide")
  if(is.null(x@mu) | is.null(x@p)) stop("Need p and mu to be set in x (use set.stats)")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")

  # check dimensions before anything
  n <- nrow(x)
  if(!missing(K)) {
    if(!is.list(K)) { 
      if(n != nrow(K) | n != ncol(K)) 
        stop("K and x dimensions don't match")
    } else {
      if(any(n != sapply(K, nrow)) | any(n != sapply(K, ncol)))
        stop("K and x dimensions don't match")
    }
  }
  if(!missing(eigenK)) {
    if(n != nrow(eigenK$vectors) | n != ncol(eigenK$vectors) | n != length(eigenK$values)) 
      stop("eigenK and x dimensions don't match")
  }

  test <- match.arg(test)
  method <- match.arg(method)

  # preparation de X 
  if(p > 0) {
    if((method == "lmm" & response == "quantitative" & test == "score") | 
       (method == "lmm" & response == "binary") |
       (method == "lm")) { # il faut ajouter les PCs à X
      X <- cbind(X, eigenK$vectors[,seq_len(p)])
      X <- trans.X(X, mean.y = mean(Y))
    } else { 
      X <- trans.X(X, eigenK$vectors[,seq_len(p)], mean(Y))
    }
  } else {
    X <- trans.X(X, mean.y = mean(Y))
  }

  # random effect
  if(method == "lmm") { 

    # if(response == "binary" & test != "score") {
    #  warning('Binary phenotype and method = "lmm" force test = "score"')
    #  test <- "score"
    # }

    if(test == "score" | response == "binary") {
      if(missing(K)) stop("For a score test and for binary traits, argument K is mandatory")
      # avec le score test on peut gérer les données manquantes dans Y
      if( any(is.na(Y)) ) {
        w <- !is.na(Y)
        X <- as.matrix(X[w,])
        Y <- Y[w]
        if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
        warning(sum(!w), 'individuals with missing phenotype are ignored.\n')
      } 
    } else {
      if(missing(eigenK)) 
        stop("For quantitative Wald and LRT tests, argument eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
    }

    if(response == "quantitative") { # score (argument K), wald ou lrt (eigen K) possibles
      if(test == "score") {
        model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
        t <- .Call(`_gaston_GWAS_lmm_score_f`, PACKAGE = "gaston", x@bed, model$Py, model$P, x@mu, beg-1, end-1)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE)
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call(`_gaston_GWAS_lmm_wald`, PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else { # test == "lrt"
        X <- cbind(X, 0) # space for the SNP
        t <- .Call(`_gaston_GWAS_lmm_lrt`, PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
      }
    } else { # response == "binary", seulement le score test, avec argument K
      if(test == "score") {
        model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
        omega <- model$BLUP_omega
        if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
        pi <- 1/(1+exp(-omega))
        t <- .Call(`_gaston_GWAS_lmm_score_f`, PACKAGE = "gaston", x@bed, Y-pi, model$P, x@mu, beg-1, end-1)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE) 
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call(`_gaston_GWAS_logitmm_wald_f`, PACKAGE = "gaston", x@bed, x@mu, Y, X, K, beg-1, end-1, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else stop("LRT test for binary trait not available")
    }
  }

  # only fixed effects
  if(method == "lm") {
    if(test != "wald") warning('Method = "lm" force test = "wald"')
    if( any(is.na(Y)) ) 
      stop("Can't handle missing data in Y")
    if(response == "quantitative") {
      t <- .Call(`_gaston_GWAS_lm_quanti`, PACKAGE = "gaston", x@bed, x@mu, Y, X, beg-1, end-1);
      t$p <- pt( abs(t$beta/t$sd), df = length(Y) - ncol(X) - 1, lower.tail=FALSE)*2
    }
    if(response == "binary") {
      X <- cbind(X,0)
      t <- .Call(`_gaston_GWAS_logit_wald_f`, PACKAGE = "gaston", x@bed, x@mu, Y, X, beg-1, end-1, tol);
      t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
    }
  }
  L <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id = x@snps$id, A1 = x@snps$A1, A2 = x@snps$A2, freqA2 = x@p)
  if(beg > 1 | end < ncol(x))  # avoid copy
    L <- L[beg:end,] 

  data.frame( c( L, t) )
}


# Check AND replaces by QR decomposition...
trans.X <- function(X, PCs = matrix(0, nrow=nrow(X), ncol=0), mean.y = 1) {
  if(any(is.na(X)))
    stop("Covariates can't be NA")

  PCs <- as.matrix(PCs) # cas où p = 1
  n.X  <- ncol(X)
  n.pc <- ncol(PCs)
  n <- n.X + n.pc

  qr.X <- qr( cbind(PCs, X) );
  if(qr.X$rank < n) {
    warning("Covariate matrix X is not full rank, removing col(s) ", paste(qr.X$pivot[ seq(qr.X$rank+1,n) ] - n.pc , collapse = ", "))
    X <- X[ , qr.X$pivot[seq(n.pc+1, qr.X$rank)] - n.pc]
    qr.X <- qr(X)
  }
  if(abs(mean.y) > 1e-4) {
    X1 <- cbind(1,X);
    qr.X1 <- qr(X1);
    if(qr.X1$rank == ncol(X1)) {
      warning("An intercept column was added to the covariate matrix X")
      X <- X1;
      qr.X <- qr.X1
    }
  }
  if( qr.X$rank == ncol(X) )
    qr.Q(qr.X)
  else
    qr.Q(qr(X))
}

