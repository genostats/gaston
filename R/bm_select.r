select.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'SNP(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[,w]
}

select.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'individual(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[w,]
}

test.snps <- function(x, condition, na.to.false = TRUE) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  miss <- is.na(w)
  if(na.to.false & sum(miss)>0)
    return(w & !miss)
  w
}

test.inds <- function(x, condition, na.to.false = TRUE) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0)
    return(w & !miss)
  w
}

which.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  which(eval(substitute(condition), x@snps, parent.frame()))
}

which.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  which(eval(substitute(condition), x@ped, parent.frame()))
}

is.autosome <- function(chr) chr %in% getOption("gaston.autosomes")
is.chr.x    <- function(chr) chr %in% getOption("gaston.chr.x")
is.chr.y    <- function(chr) chr %in% getOption("gaston.chr.y")
is.chr.mt   <- function(chr) chr %in% getOption("gaston.chr.mt")


