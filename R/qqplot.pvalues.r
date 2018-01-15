
qqplot.pvalues <- function(p, pch = ".", cex = 2, xlab = expression(paste("expected ", -log[10](p))), 
                           ylab = expression(paste("observed ", -log[10](p))), main = "QQ plot of p-values", 
                           col.abline = "red", CB = TRUE, col.CB = gray(.9), CB.level = 0.95, ...) {
  args <- list(...)
  args$pch <- pch
  args$xlab <- xlab
  args$ylab <- ylab
  args$main <- main

  n <- length(p)
  expected <- -log10( (1:n)/(n+1) )
  observed <- -log10( sort(p) )

  args$type <- "n"
  args$x <- range(expected)
  args$y <- range(observed)
  do.call( plot, args )

  # confidence interval
  k <- n * 10**seq(-log10(n), 0, length = 200)
  hi <- -log10(qbeta( 0.5-CB.level/2, k, n+1-k))
  lo <- -log10(qbeta( 0.5+CB.level/2, k, n+1-k))
  polygon( -log10(c(k, rev(k))/(n+1)), c(lo, rev(hi)), col = col.CB, border= col.CB )

  # diag line with real limits
  segments(0, 0, -log10(1/n), -log10(1/n), col = col.abline) 

  # the points
  args$x <- expected
  args$y <- observed
  args$type <- "p"
  do.call( points, args )
}

