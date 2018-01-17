
# seek duplicated by chr pos only for this version
SNP.duplicated <- function(x, by = "chr:pos") {
  if(class(x) == "bed.matrix") {
    x <- x@snps
  }
  if(class(x) != "data.frame")
    stop("x be a bed matrix or a data frame")
 
  .Call('gg_which_duplicated_chr_pos', PACKAGE = "gaston", x$chr, x$pos)
}
