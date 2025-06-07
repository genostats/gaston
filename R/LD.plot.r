#' Plot Linkage Disequilibrium
#' 
#' @description  Pretty plot of a Linkage Disequilibrium (LD) matrix
#' 
#' @param LD A symmetric LD matrix (such as produced by \code{LD}
#' @param snp.positions A vector of SNP positions
#' @param max.dist Maximal distance above which the LD is not plotted
#' @param depth Maximal number of neighbouring SNPs for which the LD is plotted
#' @param graphical.par A list of graphical parameters for function \code{par}
#' @param cex.ld The magnification to be used for LD values (if missing, an ad-hoc value is computed)
#' @param cex.snp The magnification to be used for SNPs ids (if missing, an ad-hoc value is computed)
#' @param polygon.par A list of parameters for function \code{polygon}
#' @param color.scheme A function to set the background color of a cell
#' @param write.snp.id \code{Logical}. If \code{TRUE}, SNP ids will be displayed above the plot
#' @param write.ld \code{NULL}, or a function which outputs the string used for displaying a LD value in a cell
#' @param draw.chr \code{Logical}. If \code{TRUE}, a chromosome with SNP positions is sketched above the plot
#' @param above.space Space above the plot (in user units = height of a cell)
#' @param below.space Space below the plot (in user units = height of a cell)
#' @param pdf.file The name of a pdf file in which to plot the LD matrix. If missing, current plot device will be used
#' @param finalize.pdf \code{Logical}. If \code{TRUE}, \code{dev.off()} will be called to finalize the pdf file
#' 
#' 
#' @details This function displays a LD plot similar to Haploview plots.
#' 
#' To add anotations to the plot, it is useful to know that each cell has width and height equal
#' to one user unit, the first cell in the upper row being centered at coordinates \code{(1.5, -0.5)}.
#' 
#' @seealso  \code{\link{LD}}
#' 
#' @keywords  Linkage Disequilibrium
#' @examples
#' 
#' # Load data
#' data(AGT)
#' x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
#' 
#' # Compute LD
#' ld.x <- LD(x, c(1,ncol(x)))
#' 
#' # Plot a tiny part of the LD matrix
#' LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )
#' 
#' # Customize the plot
#' LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20], 
#'          graphical.par = list(cex = 1.3, bg = "gray"), 
#'          polygon.par = list(border = NA), write.ld = NULL )
#' \dontrun{
#' # Plotting the whole matrix in X11 display is very long (lots of polygons)
#' # but it is ok with a pdf file
#' # (please uncomment to run)
#' #LD.plot(ld.x, snp.positions = x@snps$pos, max.dist = 50e3, write.ld = NULL, pdf.file = "LDAGT.pdf")
#' }
#' 
#' @export LD.plot
LD.plot <- function(LD, snp.positions, max.dist = Inf, depth = nrow(LD), graphical.par = list(mar = c(0,0,0,0)),
                    cex.ld, cex.snp, polygon.par = list(border = "white"),
                    color.scheme = function(ld) rgb(1,1-abs(ld),1-abs(ld)),
                    write.snp.id = TRUE, write.ld = function(ld) sprintf("%.2f", ld),
                    draw.chr = TRUE,
                    above.space = 1 + 2*write.snp.id + draw.chr, below.space = 1, 
                    pdf.file, finalize.pdf = TRUE) {

  n <- nrow(LD)
  positions <- if(missing(snp.positions)) rep(0,n) else snp.positions
  # dry run to get graph depth
  graph.depth <- 0
  for(i in seq(1,n-1)) 
    for(j in seq(i+1, min(n,i+depth))) {
      if(positions[j] - positions[i] > max.dist) next;
      graph.depth <- max(graph.depth, (j-i)/2)
    }


  if(!missing(pdf.file)) 
    pdf(pdf.file, width = n/2, height = (graph.depth + above.space + below.space)/2)

  do.call(par, graphical.par) 
  write.ld.val <- ifelse(is.null(write.ld), FALSE, TRUE)
  plot( 0,0, xlim=c(0.5,n+0.5), ylim=c(-graph.depth-below.space,above.space), type="n", 
        asp = 1, xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  ld <- ""

  if(missing(cex.ld) & write.ld.val) 
    cex.ld = 0.7*par("cex") / max(strwidth(sapply(LD, write.ld)))

  for(i in seq(1,n-1)) 
    for(j in seq(i+1, min(n,i+depth))) {
      if(positions[j] - positions[i] > max.dist) next;
      if(write.ld.val) ld <- write.ld(LD[i,j])
      snp.sq(i, j, ld, write.ld.val, color.scheme(LD[i,j]), polygon.par, cex.ld)
    }
  if(write.snp.id) {
    rs.h <- strheight(rownames(LD))
    if(missing(cex.snp)) cex.snp <- 0.25*par("cex")/max(rs.h)
    text( 1:n, 0, rownames(LD), srt = 90, cex = cex.snp, adj = c(0,0.5))
  }
  if(!missing(snp.positions) & draw.chr) {
    if(write.snp.id) 
      a <- max( strwidth(rownames(LD), cex = cex.snp) )
    else
      a <- 0
    pos <- 1.5 + (n-2)*(snp.positions - snp.positions[1])/diff(range(snp.positions))
    segments( 1:n, a+0.25,  pos, a + 1.5 )
    rect(1.5, a+1.5, n-0.5, a+1.75)
    segments( pos, a+1.5,  pos, a + 1.75 )
  }
  if(!missing(pdf.file) & finalize.pdf) dev.off()
}

# SNP square
snp.sq <- function(i, j, ld, write.ld.val, color, polygon.par, cex.ld) {
  if(i == j) return(); # pas de plot d'un SNP avec lui même
  if(j < i) { tmp <- j; j <- i; i <- tmp } # i < j
  d <- (j-i)/2
  cx <- i+d
  cy <- -d
  do.call(polygon, c(list(x = cx+c(-1,0,1,0)/2, y = cy + c(0,1,0,-1)/2, col = color), polygon.par))
  if(write.ld.val) text(cx, cy, ld, cex = cex.ld)
}



