
SNP.match <- function(x, table, by = "chr:pos:alleles") {
  # préparation de x = data frame avec les bonnes colonnes
  b <- strsplit(by,':')[[1]]
  if ('alleles' %in% b) 
    b <- c(b[b!='alleles'], 'A1', 'A2')

  if (is(x, "bed.matrix")) {
    x <- x@snps
  }

  x <- x[, b, drop = FALSE]
  # c'est prêt

  if (is(table, "bed.matrix")) {
    table <- table@snps
  }
  if (!inherits(x, "data.frame") || !inherits(table, "data.frame"))
    stop("x and table should be bed matrices or data frames")
  
  .Call('gg_SNPmatch', PACKAGE = "gaston", x, table)
}
