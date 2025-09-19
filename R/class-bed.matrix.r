#' Class \code{"bed.matrix"}
#'
#' @name bed.matrix
#' @aliases bed.matrix-class
#' @docType class
#' 
#' @description  S4 class for SNP genotype matrices
#' 
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("bed.matrix", ...)}.
#'
#' @slot ped 
#' \code{data.frame} containing information for each individual: \code{famid} = Family ID,
#' \code{id} = Individual ID, \code{father} = Father ID, \code{mother} = Mother ID, \code{sex} = Sex and \code{pheno} = Phenotype.
#' Can also contain individuals statistic, for example: \code{N0}, \code{N1} and \code{N2} = Number of genotypes equal to 0, 1 and 2 respectively,
#' \code{NAs} = Number of missing genotypes, \code{callrate} = Individual callrate.
#' @slot snps 
#' \code{data.frame} containing information for each SNP: \code{chr} = Chromosome, \code{id} = SNP ID,
#' \code{dist} = Genetic Distance, \code{pos} = Physical position, \code{A1} = Reference Allele, \code{A2} = Alternative Allele.
#' Can also contain SNPs statistic, for example: \code{N0}, \code{N1} and \code{N2} = Number of genotypes equal to 0, 1 and 2 repectively,
#' \code{NAs} = Number of missing genotypes, \code{callrate} = SNP callrate, \code{maf} = Minor allele frequency), \code{hz} = heterozygosity
#' @slot bed 
#' \code{externalptr} (pointing to the genotype matrix). 
#' @slot p
#' \code{vector} or \code{NULL} for allelic frequencies (allèle \code{A2}).
#' @slot mu 
#' \code{vector} or \code{NULL} for genotype means (usually \code{mu = 2*p}).
#' @slot sigma
#' \code{vector} or \code{NULL} for genotypic standard deviation
#' @slot standardize_p
#' \code{logical}. If \code{TRUE}, the genotype matrix is standardized using means \code{2*p}
#' and genotypic standard deviation \code{sqrt(2*p*(1-p))}
#' @slot standardize_mu_sigma
#' \code{logical}. If \code{TRUE}, the genotype matrix is standardize using means
#' \code{mu} and genotypic standard deviation \code{sigma}.
#' For more details please check the vignette.
#' @section Methods:
#' \describe{
#' \item{as.matrix}{\code{signature(x = "bed.matrix")}:
#' \cr Convert a \code{bed.matrix} into a \code{matrix}.
#' The resulting matrix can be huge, use this method only for a small bed.matrix! }
#' \item{standardize}{\code{signature(x = "bed.matrix")}:
#' \cr Get the standardize parameter of \code{bed.matrix}. Can be "none", "p" or "mu_sigma". }
#' \item{standardize<-}{\code{signature(x = "bed.matrix")}:
#' \cr Set the standardize parameter of a \code{bed.matrix}. }
#' \item{dim}{\code{signature(x = "bed.matrix")}:
#' \cr Get the number of individuals (rows) and the number of SNPs (columns). }
#' \item{head}{\code{signature(x = "bed.matrix")}:
#' \cr Print the head of the genotype matrix of a \code{bed.matrix} object. }
#' \item{mu}{\code{signature(x = "bed.matrix")}:
#' \cr Get the \code{mu} slot of a \code{bed.matrix}. }
#' \item{mu<-}{\code{signature(x = "bed.matrix")}:
#' \cr Set the \code{mu} slot of a \code{bed.matrix}. }
#' \item{p}{\code{signature(x = "bed.matrix")}:
#' \cr Get the \code{p} slot of a \code{bed.matrix}. }
#' \item{p<-}{\code{signature(x = "bed.matrix")}:
#' \cr Set the \code{p} slot of a \code{bed.matrix}. }
#' \item{show}{\code{signature(object = "bed.matrix")}:
#' \cr Displays basic information about a \code{bed.matrix}. }
#' \item{sigma}{\code{signature(x = "bed.matrix")}:
#' \cr Get the \code{sigma} slot of a \code{bed.matrix}. }
#' \item{sigma<-}{\code{signature(x = "bed.matrix")}:
#' \cr Set the \code{sigma} slot of a \code{bed.matrix}. }
#' \item{cbind}{\code{signature(... = "bed.matrix")}:
#' \cr Combine a sequence of \code{bed.matrix} by columns. }
#' \item{rbind}{\code{signature(... = "bed.matrix")}:
#' \cr Combine a sequence of \code{bed.matrix} by rows. }
#' }
#' 
#' @seealso  \code{\link{read.bed.matrix}}, \code{\link{write.bed.matrix}},
#' \code{\link{set.stats}}, \code{\link{select.snps}}, \code{\link{select.inds}}, \code{\link{GRM}}
#' 
#' @keywords classes
#' @examples
#' 
#' showClass("bed.matrix")
#' 
#' # Conversion example
#' data(LCT)
#' x1 <- as(LCT.gen, "bed.matrix")
#' x1
#' head(x1@ped)
#' head(x1@snps)
#' 
#' # the function as.bed.matrix is an alternative
#' x2 <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
#' x2
#' head(x2@ped)
#' head(x2@snps)
#' 

setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))

pednames <- c("famid", "id", "father", "mother", "sex", "pheno")
snpnames <- c("chr", "id", "dist", "pos", "A1", "A2")

snpstatnames0 <- c("N0", "N1", "N2", "NAs", "callrate", "maf", "hz", "N0.f", "N1.f", "N2.f", "NAs.f" )
snpstatnames <- c(snpstatnames0, "hwe")
pedstatnames <- c("N0", "N1", "N2", "NAs", "N0.x", "N1.x", "N2.x", "NAs.x", 
                  "N0.y", "N1.y", "N2.y", "NAs.y", "N0.mt", "N1.mt", "N2.mt", "NAs.mt", 
                  "callrate", "hz", "callrate.x", "hz.x", "callrate.y", "hz.y", "callrate.mt", "hz.mt")

is.null.df <- function(x) is.data.frame(x) & nrow(x) == 0 & ncol(x) == 0

## Class bed.matrix
#' @exportClass bed.matrix
setClass("bed.matrix", representation( 
                   ped = 'data.frame',
                   snps = 'data.frame', 
                   bed = 'externalptr',
                   p = 'numericOrNULL',
                   mu = 'numericOrNULL',
                   sigma = 'numericOrNULL',
                   standardize_p = 'logical',
                   standardize_mu_sigma = 'logical' ))

setValidity('bed.matrix',
           function(object) {
             errors <- character()
             if ( object@standardize_p & object@standardize_mu_sigma ) 
                errors <- c(errors, "Only one center scale parameter can be TRUE.")
             if ( object@standardize_p & is.null(object@p) ) 
                errors <- c(errors, "If 'standardize_p' is TRUE, 'p' must be defined.")
             if ( object@standardize_mu_sigma & ( is.null(object@mu) | is.null(object@sigma) ) ) 
                errors <- c(errors, "If 'standardize_mu_sigma' is TRUE, 'mu' and 'sigma' must be defined.")
             if ( !is.null(object@p) & length(object@p) != ncol(object) ) 
                errors <- c(errors, "The length of 'p' must be equal to the number of markers.")
             if ( !is.null(object@mu) & length(object@mu) != ncol(object) ) 
                errors <- c(errors, "The length of 'mu' must be equal to the number of markers.")
             if ( !is.null(object@sigma) & length(object@sigma) != ncol(object) ) 
                errors <- c(errors, "The length of 'sigma' must be equal to the number of markers.")
             if(.Call(`_gaston_isnullptr`,  PACKAGE = 'gaston', object@bed))
                errors <- c(errors, 'The externalptr is broken')
             if ( length(errors)==0 ) return(TRUE) else return(errors)
           } );


setAs("bed.matrix", "matrix",
  function(from) {
    validObject(from)
    to <- if(from@standardize_p) 
      .Call(`_gaston_m4_as_scaled_matrix_p`, PACKAGE = 'gaston', from@bed, from@p)
    else if(from@standardize_mu_sigma)
      .Call(`_gaston_m4_as_scaled_matrix_mu_sigma`, PACKAGE = 'gaston', from@bed, from@mu, from@sigma)
    else
      .Call(`_gaston_m4_as012`, PACKAGE = 'gaston', from@bed)
    colnames(to) <- from@snps$id
    rownames(to) <- if(any(duplicated(from@ped$id))) paste(from@ped$fam, from@ped$id, sep=":")
                    else from@ped$id
    to
  } );

setGeneric('as.matrix')

#' @exportMethod as.matrix
setMethod("as.matrix", signature="bed.matrix", definition = function(x) as(x,"matrix") )

setAs("matrix", "bed.matrix", 
  function(from) {
    bed <- .Call(`_gaston_as_matrix4`, PACKAGE = 'gaston', from)

    ped <- if(is.null(rownames(from))) 
             structure(list(), row.names = c(NA, -nrow(from)), class = "data.frame") # empty data frame with right number of lines
           else 
             data.frame(famid = rownames(from), id = rownames(from), father = 0, mother = 0, sex = 0, pheno = 0, stringsAsFactors = FALSE)

    snp <- if(is.null(colnames(from)))
             structure(list(), row.names = c(NA, -ncol(from)), class = "data.frame") #idem
           else 
             data.frame(chr = NA, id = colnames(from), dist = NA, pos = NA, A1 = NA, A2 = NA, stringsAsFactors = FALSE)

    x <- new("bed.matrix", bed = bed, snps = snp, ped = ped, p = NULL, mu = NULL,
             sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
    if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
    x
  } );

setGeneric('dim')

#' @exportMethod dim
setMethod("dim", signature = "bed.matrix", 
    function(x) c(.Call(`_gaston_ninds`, PACKAGE = 'gaston', x@bed), .Call(`_gaston_nsnps`, PACKAGE = 'gaston', x@bed)))

setGeneric('head')

#' @exportMethod head
setMethod( 'head', signature(x='bed.matrix'), function(x, nrow=10, ncol=10) print( as.matrix( x[1:min( nrow, nrow(x) ),1:min( ncol, ncol(x) )] ) ) )

setMethod(show, signature("bed.matrix"), 
  function(object) {
    if(.Call(`_gaston_isnullptr`,  PACKAGE = 'gaston', object@bed))
      cat("A bed.matrix with a broken externalptr!\nHint: don't save/load bed.matrices with other functions than write.bed.matrix and read.bed.matrix\n")
    else {
      cat('A bed.matrix with ', nrow(object), ' individuals and ', ncol(object), ' markers.\n', sep='')
      if(anyDuplicated(object@snps$id)) cat("There are some duplicated SNP id's\n")
      if(all(snpstatnames0 %in% names(object@snps))) {
        cat("snps stats are set\n");
        a <- sum(object@snps$NAs == nrow(object))
        if(a > 1)  cat("  There are ", a, " SNPs with a callrate equal to zero\n");
        if(a == 1) cat("  There is one SNP with a callrate equal to zero\n");
        a <- sum(object@snps$maf == 0, na.rm = TRUE)
        if(a > 1)  cat("  There are ", a, " monomorphic SNPs\n");
        if(a == 1) cat("  There is one monomorphic SNP\n");
      } else 
        cat("snps stats are not set (or incomplete)\n")
      # ped stats et autres
      if(anyDuplicated(object@ped[, c("famid", "id")])) 
        cat("There are some duplicated individual id's\n")
      if(all(pedstatnames %in% names(object@ped))) {
        cat("ped stats are set\n");
        a <- sum(object@ped$NAs == ncol(object))
        if(a > 1)  cat("  There are ", a, " individuals with a callrate equal to zero\n");
        if(a == 1) cat("  There is one individual with a callrate equal to zero\n");
      } else cat("ped stats are not set (or incomplete)\n")
    }
  } 
)

