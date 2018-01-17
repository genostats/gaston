
# seek duplicated by chr pos only for this version
SNP.duplicated <- function(x, by = "chr:pos") {
  if(class(x) == "bed.matrix") {
    x <- x@snps
  }
  if(class(x) != "data.frame")
    stop("x be a bed matrix or a data frame")

  if(by == "id") {
    if(is.null(x$id)) stop("Missing id in x")
    return(.Call('gg_which_duplicated_id', PACKAGE = "gaston", x$id))
  }
  if(by == "id:chr:pos") {
    if(is.null(x$id)) stop("Missing id in x")
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    return(.Call('gg_which_duplicated_id_chr_pos', PACKAGE = "gaston", x$id, x$chr, x$pos))
  }

  if(by == "id:chr:pos:alleles") {
    if(is.null(x$id)) stop("Missing id in x")
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    if(is.null(x$A1) || is.null(x$A2)) stop("Missing alleles in x")
    return(.Call('gg_which_duplicated_id_chr_pos_alleles', PACKAGE = "gaston", x$id, x$chr, x$pos, x$A1, x$A2))
  }

  if(by == "chr:pos") {
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    return(.Call('gg_which_duplicated_chr_pos', PACKAGE = "gaston", x$chr, x$pos))
  }

  if(by == "chr:pos:alleles") {
    if(is.null(x$chr)) stop("Missing chr in x")
    if(is.null(x$pos)) stop("Missing pos in x")
    if(is.null(x$A1) || is.null(x$A2)) stop("Missing alleles in x")
    return(.Call('gg_which_duplicated_chr_pos_alleles', PACKAGE = "gaston", x$chr, x$pos, x$A1, x$A2))
  }
  stop("'by' should be one of 'id', 'id:chr:pos', 'id:chr:pos:alleles', 'chr:pos', 'chr:pos:alleles'")
}
