#' Manhattan plot
#' 
#' @description  Draws a Manhattan plot
#' 
#' @param x  A data.frame with columns named \code{chr}, \code{pos} and \code{p}.
#' @param bty  Type of box to draw about the plot. Default is to draw none.
#' @param thinning  \code{Logical}. If \code{TRUE}, not all points are displayed.
#' @param chrom.col  Alternating colors for chromosomes.
#' @param ...  Graphical parameters to be passed to \code{plot}.
#' 
#' 
#' @details If there is only one chromosome value in \code{x$chr}, the x-axis will be labeled with the SNP
#' position. In the general case, the x-axis is labeled with the chromosome name and the color
#' of the points alternates between the colors in \code{chrom.col}.
#' 
#' The default value \code{bty = "n"} should give the best result for GWAS Manhattan plots.
#' See \code{\link{par}} for other possible values of \code{bty} and their meaning.
#' 
#' The thinning procedure suppress some points to avoid generating too heavy graphs. The user
#' should check that setting \code{thinning = FALSE} does not change the final aspect of the
#' plot.
#' 
#' 
#' 
#' @return An invisible copy of \code{x} is returned, in which a column \code{coord} has been added
#' if there is more than one chromosome value in \code{x$chr}. This column contains the x-coordinates of each
#' SNP on the plot, and should prove helpful to annotate it.
#' @seealso  \code{\link{association.test}}, \code{\link{qqplot.pvalues}},
#' \code{\link{par}}, \code{\link{plot.default}}, \code{\link{points.default}}
#' 
#' 
#' @export manhattan
manhattan <- function(x, bty="n", chrom.col = c("black", "gray50"), thinning = TRUE, ... ) {
  if(!is.data.frame(x)) 
    stop("x should be a data.frame")
  if(!all( c("chr", "pos", "p") %in% names(x) )) 
    stop("x should have 'chr', 'pos' and 'p' components")


  tp0 <- (x$p == 0)
  if(any(tp0, na.rm = TRUE)) warning("There are ", sum(tp0, na.rm = TRUE), " null p-values that won't appear on the plot")

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
