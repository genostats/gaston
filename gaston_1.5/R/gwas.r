association.test <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), 
                             method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                             test = c("score", "wald", "lrt"), 
                             K, eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, multithreaded = FALSE, ...) {

  if(beg < 1 || end > ncol(x)) stop("range too wide")
  if(is.null(x@mu)) stop("Need mu to be set in x")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")

  response <- match.arg(response)
  test <- match.arg(test)

  # random effect
  if(match.arg(method) == "lmm") { # il n'y a que ça pour le moment

    if(response == "binary" & test != "score") {
      warning('Binary phenotype and method = "lmm" force test = "score"')
      test <- "score"
    }

    if(test == "score") {
      if(missing(K)) stop("For a score test, argument K is mandatory")
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
        stop("For Wald and LRT tests, argument eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
    }

    if(response == "quantitative") { # score (argument K), wald ou lrt (eigen K) possibles
      if(test == "score") {
        model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
        t <- .Call("gg_GWAS_lmm_score", PACKAGE = 'gaston', x@bed, model$Py, model$P, x@mu, beg-1, end-1)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE)
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        if(ncol(x) > 2000 & multithreaded) 
          t <- .Call("gg_GWAS_lmm_wald_mt", PACKAGE = 'gaston', x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        else
          t <- .Call("gg_GWAS_lmm_wald", PACKAGE = 'gaston', x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else { # test == "lrt"
        X <- cbind(X, 0) # space for the SNP
        if(ncol(x) > 2000 & multithreaded) 
          t <- .Call("gg_GWAS_lmm_lrt_mt", PACKAGE = 'gaston', x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        else
          t <- .Call("gg_GWAS_lmm_lrt", PACKAGE = 'gaston', x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
        t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
      }
    } else { # response == "binary", seulement le score test, avec argument K
      model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
      omega <- model$BLUP_omega
      if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
      pi <- 1/(1+exp(-omega))
      t <- .Call("gg_GWAS_lmm_score", PACKAGE = 'gaston', x@bed, Y-pi, model$P, x@mu, beg-1, end-1)
      t$p <- pchisq( t$score, df = 1, lower.tail=FALSE)
    }
  }

  # only fixed effects
  if(match.arg(method) == "lm") {
    if(response == "quantitative") {
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      if(p > 0)
        Q <- qr.Q(qr(cbind(eigenK$vectors[,seq_len(p)], X)))
      else
        Q <- qr.Q(qr(X))
      t <- .Call("gg_GWAS_lm_quanti", PACKAGE = 'gaston', x@bed, x@mu, Y, Q, beg-1, end-1);
      t$p <- pt( abs(t$beta/t$sd), df = length(Y) - ncol(Q) - 1, lower.tail=FALSE)*2
    }
    if(response == "binary") {
      stop("logistic regression not yet available\n")
    }
  }


  return(data.frame(t))
}
