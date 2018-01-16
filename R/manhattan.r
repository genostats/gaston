manhattan <- function(x, bty="n", thinning = TRUE, thinning.step = 1e-5, ... ) {
  if(!is.data.frame(x)) 
    stop("x should be a data.frame")
  if(!all( c("chr", "pos", "p") %in% names(x) )) 
    stop("x should have 'chr', 'pos' and 'p' components")

  chr <- as.factor(x$chr)
  if(nrow(x) < 1e5) thinning <- FALSE # disable thinning for 'small' sets
  # create coordinate column
  if(nlevels(chr) > 1) {
    x$coord <- 0
    M <- 0
    tic <- numeric( nlevels(chr) )
    for(i in 1:nlevels(chr)) {
      w <- (chr == levels(chr)[i])
      pos.c <- x$pos[w]
      x$coord[w] <- M + pos.c
      mx <- max( pos.c )
      tic[i] <- M + mx/2
      M <- M + mx
    }
    x$coord <- x$coord/M
    tic <- tic/M
  }

  args <- list(...)

  if(is.null(args$xlab))
    args$xlab <- if(nlevels(chr) > 1) "Chromosome" else "Position"
  if(is.null(args$ylab))
    args$ylab <- expression(-log[10](p))
  if(is.null(args$col))
    args$col <- ifelse( as.integer(chr) %% 2, "black", "gray50" )
  else
    args$col <- rep_len(args$col, nrow(x))

  args$bty  <- bty
  args$xaxt <- "n" 
    
  args$x <- if(nlevels(chr) > 1) x$coord else x$pos
 
  if(thinning) {
    lp <- -log10(x$p)
    o <- order( lp )
    lp <- lp[o]
    w <- gaston:::logp.thinning(lp, thinning.step)
    args$y <- lp[w]
    oo <- o[w]
    args$x <- args$x[oo]
    args$col  <- args$col[oo]
    if(length(args$pch) > 1) args$pch <- rep_len(args$pch, nrow(x))[oo]
    if(length(args$cex) > 1) args$cex <- rep_len(args$cex, nrow(x))[oo]
    if(length(args$lwd) > 1) args$lwd <- rep_len(args$lwd, nrow(x))[oo]
    if(length(args$bg)  > 1) args$bg  <- rep_len(args$bg,  nrow(x))[oo]
  } else {
    args$y <- -log10(x$p) 
  }
  do.call(plot, args)
 
  # draw axis unless xaxt = "n" given by user... 
  ax <- list(...)$xaxt
  ax <- if(is.null(ax)) "s" else ax 
  if(nlevels(chr) > 1) 
    axis(1, tic, levels(chr), xaxt = ax)
  else
    axis(1, xaxt = ax) 
  invisible(x);
}
