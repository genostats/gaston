association.test <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), 
                             method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                             test = c("score", "wald", "lrt"), 
                             K, eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, ...) {

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

  response <- match.arg(response)
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
        t <- .Call("gg_GWAS_lmm_score_f", PACKAGE = "gaston", x@bed, model$Py, model$P, x@mu, beg-1, end-1)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE)
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("gg_GWAS_lmm_wald", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else { # test == "lrt"
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("gg_GWAS_lmm_lrt", PACKAGE = "gaston", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
      }
    } else { # response == "binary", seulement le score test, avec argument K
      if(test == "score") {
        model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
        omega <- model$BLUP_omega
        if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
        pi <- 1/(1+exp(-omega))
        t <- .Call("gg_GWAS_lmm_score_f", PACKAGE = "gaston", x@bed, Y-pi, model$P, x@mu, beg-1, end-1)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE) 
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("gg_GWAS_logitmm_wald_f", PACKAGE = "gaston", x@bed, x@mu, Y, X, K, beg-1, end-1, tol)
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
      t <- .Call("gg_GWAS_lm_quanti", PACKAGE = "gaston", x@bed, x@mu, Y, X, beg-1, end-1);
      t$p <- pt( abs(t$beta/t$sd), df = length(Y) - ncol(X) - 1, lower.tail=FALSE)*2
    }
    if(response == "binary") {
      X <- cbind(X,0)
      t <- .Call("gg_GWAS_logit_wald_f", PACKAGE = "gaston", x@bed, x@mu, Y, X, beg-1, end-1, tol);
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
  if(mean.y > 1e-4) {
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

