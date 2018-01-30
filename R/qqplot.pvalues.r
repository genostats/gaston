qqplot.pvalues <- function(p, col.abline = "red", CB = TRUE, col.CB = "gray80", CB.level = 0.95,
                           thinning = TRUE, ...) {
  if(is.list(p)) { # ok pour les data frame aussi
    if(is.null(p$p))
      stop("No p-values were found")
    p <- p$p
  }
  p <- p[ !is.na(p) ] # on supprime ces valeurs silencieusement...
  
  if(any(p>1) | any(p<0))
    stop("p-values should be in [0,1]\n")

  w <- (p == 0)
  if(any(w)) { # mais ça on avertit
    warning("There are ", sum(w), " zero p-values that won't be displayed")
    p <- p[!w]
  }

  args <- list(...)
  if(is.null(args$xlab))
    args$xlab <- expression(paste("expected ", -log[10](p)))
  if(is.null(args$ylab))
    args$ylab <- expression(paste("observed ", -log[10](p)))
  if(is.null(args$main))
    args$main <- "QQ plot of p-values"

  n <- length(p)
  args$type <- "n"
  args$x <- -log10( c(1,n) / (n+1) )
  args$y <- -log10( range(p) )
  do.call( plot, args )

  # confidence interval
  k <- n * 10**seq(-log10(n), 0, length = 200)
  hi <- -log10(qbeta( 0.5-CB.level/2, k, n+1-k))
  lo <- -log10(qbeta( 0.5+CB.level/2, k, n+1-k))
  polygon( -log10(c(k, rev(k))/(n+1)), c(lo, rev(hi)), col = col.CB, border= col.CB )

  # diag line with real limits
  segments(0, 0, -log10(1/n), -log10(1/n), col = col.abline) 

  # the points
  expected <- -log10( (n:1)/(n+1) )
  observed <- sort(-log10( p ))

  if(thinning) {
    w <- manhattan.thinning(expected, observed, 10000, 10000)
    args$x <- expected[w]
    args$y <- observed[w] 
  } else {
    args$x <- expected
    args$y <- observed
  }
  args$type <- "p"
  do.call( points, args )
}

