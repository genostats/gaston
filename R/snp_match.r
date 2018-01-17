
SNP.match <- function(x, table) {
  if(class(x) == "bed.matrix") {
    x <- x@snps
  }
  if(class(table) == "bed.matrix") {
    table <- table@snps
  }
  if(class(x) != "data.frame" | class(table) != "data.frame")
    stop("x and table should be bed matrices or data frames")
  
  .Call('gg_SNPmatch', PACKAGE = "gaston", x, table)
}
