\name{LD.plot}
\alias{LD.plot}

\title{ Plot Linkage Disequilibrium }

\description{ Pretty plot of a Linkage Disequilibrium (LD) matrix }
\usage{ LD.plot(LD, snp.positions, max.dist = Inf, depth = nrow(LD), 
                graphical.par = list(mar = c(0,0,0,0)), cex.ld, cex.snp,
                polygon.par = list(border = "white"), 
                color.scheme = function(ld) rgb(1,1-abs(ld),1-abs(ld)),
                write.snp.id = TRUE, write.ld = function(ld) sprintf("\%.2f", ld),
                draw.chr = TRUE, above.space = 1 + 2*write.snp.id + draw.chr, 
                below.space = 1, pdf.file, finalize.pdf = TRUE) 
}
\arguments{
  \item{LD}{A symmetric LD matrix (such as produced by \code{LD}}
  \item{snp.positions}{A vector of SNP positions}
  \item{max.dist}{Maximal distance above which the LD is not plotted}
  \item{depth}{Maximal number of neighbouring SNPs for which the LD is plotted}  
  \item{graphical.par}{A list of graphical parameters for function \code{par}}  
  \item{cex.ld}{The magnification to be used for LD values (if missing, an ad-hoc value is computed)}  
  \item{cex.snp}{The magnification to be used for SNPs ids (if missing, an ad-hoc value is computed)}  
  \item{polygon.par}{A list of parameters for function \code{polygon}}  
  \item{color.scheme}{A function to set the background color of a cell}  
  \item{write.snp.id}{\code{Logical}. If \code{TRUE}, SNP ids will be displayed above the plot}  
  \item{write.ld}{\code{NULL}, or a function which outputs the string used for displaying a LD value in a cell}  
  \item{draw.chr}{\code{Logical}. If \code{TRUE}, a chromosome with SNP positions is sketched above the plot }  
  \item{above.space}{Space above the plot (in user units = height of a cell)}  
  \item{below.space}{Space below the plot (in user units = height of a cell)}  
  \item{pdf.file}{The name of a pdf file in which to plot the LD matrix. If missing, current plot device will be used}  
  \item{finalize.pdf}{\code{Logical}. If \code{TRUE}, \code{dev.off()} will be called to finalize the pdf file} 
}
\details{This function displays a LD plot similar to Haploview plots. 

To add anotations to the plot, it is useful to know that each cell has width and height equal 
to one user unit, the first cell in the upper row being centered at coordinates \code{(1.5, -0.5)}.
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{LD}} }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)

# Compute LD
ld.x <- LD(x, c(1,ncol(x)))

# Plot a tiny part of the LD matrix
LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20] )

# Customize the plot
LD.plot( ld.x[1:20,1:20], snp.positions = x@snps$pos[1:20], 
         graphical.par = list(cex = 1.3, bg = "gray"), 
         polygon.par = list(border = NA), write.ld = NULL )
\dontrun{
# Plotting the whole matrix in X11 display is very long (lots of polygons)
# but it is ok with a pdf file
# (please uncomment to run)
#LD.plot(ld.x, snp.positions = x@snps$pos, max.dist = 50e3, write.ld = NULL, pdf.file = "LDAGT.pdf")
}
}

\keyword{ Linkage Disequilibrium }
