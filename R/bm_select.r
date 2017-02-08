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

is.autosome <- function(chr) chr %in% getOption("gaston.autosomes")
is.chr.x    <- function(chr) chr %in% getOption("gaston.chr.x")
is.chr.y    <- function(chr) chr %in% getOption("gaston.chr.y")
is.chr.mt   <- function(chr) chr %in% getOption("gaston.chr.mt")
