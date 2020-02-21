GRM <- function(x, which.snps, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  if(!x@standardize_mu_sigma & !x@standardize_p) {
    if(!is.null(x@p)) x@standardize_p <- TRUE
    else stop("Can't standardize x for GRM computation (use set.stat)\n")
  }

  if(x@standardize_mu_sigma) {
    w <- ifelse(x@sigma == 0, 0, 1/x@sigma/sqrt(sum(which.snps)-1))   ### BEWARE q-1 !!!
    K <- .Call('gg_Kinship_w', PACKAGE = "gaston", x@bed, x@mu[which.snps], w[which.snps], which.snps, chunk) 
  } else { 
    K <- .Call('gg_Kinship_pw', PACKAGE = "gaston", x@bed, x@p[which.snps], which.snps, FALSE, chunk)
  }

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0)
      rownames(K) <- colnames(K) <- x@ped$id
    else {
      nn <- paste(x@ped$famid, x@ped$id, sep = ":")
      if(anyDuplicated(nn) == 0)
        rownames(K) <- colnames(K) <- nn
    }
  }

  K
}

DM <- function(x, which.snps, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  if(x@standardize_mu_sigma)
    warning("For Dominance Matrix, p standardization is used\n")
  if(is.null(x@p))
    stop("Can't standardize x for DM computation (use set.stat)\n")

  K <- .Call('gg_Kinship_pw', PACKAGE = "gaston", x@bed, x@p[which.snps], which.snps, TRUE, chunk)

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K) <- colnames(K) <- x@ped$id
    else {
      nn <- paste(x@ped$famid, x@ped$id, sep = ":")
      if(anyDuplicated(nn) == 0)
        rownames(K) <- colnames(K) <- nn
    }
  }

  K
}

reshape.GRM <- function(K, include = c(-Inf, +Inf), exclude) {
  diag(K) <- NA
  if(missing(exclude))
    w <- which(include[1] < K & K < include[2])
  else 
    w <- which(include[1] < K & K < include[2] & (K < exclude[1] | K > exclude[2]))
  I <- row(K)[w]
  J <- col(K)[w]
  R <- K[w]
  ww <- (I < J)
  i <- I[ww];
  j <- J[ww];
  data.frame(i = i, j = j, id_i = rownames(K)[i], id_j = colnames(K)[j], k = R[ww])
}

