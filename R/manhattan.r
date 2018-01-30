manhattan <- function(x, bty="n", chrom.col = c("black", "gray50"), thinning = TRUE, ... ) {
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
    if(nlevels(chr) > 1)
      args$col <- chrom.col[ 1 + as.integer(chr) %% length(chrom.col) ]
    else 
      args$col <- chrom.col[1]
  else if(length(args$col) > 1) 
      args$col <- rep_len(args$col, nrow(x))

  args$bty  <- bty
  args$xaxt <- "n" 
    
  args$x <- if(nlevels(chr) > 1) x$coord else x$pos
  args$y <- -log10(x$p) 
 
  if(thinning) {
    # ces valeurs ont l'air à peu près adaptées à n'importe quel manhattan plot, 
    # pas la peine d'enquiquiner l'utilisateur avec ces paramètres (en plus il faudrait les
    # documenter)
    w <- manhattan.thinning(args$x, args$y, 1000, 1000) 
    args$y <- args$y[w]
    args$x <- args$x[w]
    if(length(args$col) > 1) args$col <- args$col[w]
    if(length(args$pch) > 1) args$pch <- rep_len(args$pch, nrow(x))[w]
    if(length(args$cex) > 1) args$cex <- rep_len(args$cex, nrow(x))[w]
    if(length(args$lwd) > 1) args$lwd <- rep_len(args$lwd, nrow(x))[w]
    if(length(args$bg)  > 1) args$bg  <- rep_len(args$bg,  nrow(x))[w]
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
