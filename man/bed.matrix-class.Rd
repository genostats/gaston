\name{bed.matrix-class}
\Rdversion{1.1}
\docType{class}
\alias{bed.matrix-class}
\alias{bed.matrix}
\alias{[,bed.matrix,logical,logical,missing-method}
\alias{[,bed.matrix,logical,missing,missing-method}
\alias{[,bed.matrix,logical,numeric,missing-method}
\alias{[,bed.matrix,missing,logical,missing-method}
\alias{[,bed.matrix,missing,numeric,missing-method}
\alias{[,bed.matrix,numeric,logical,missing-method}
\alias{[,bed.matrix,numeric,missing,missing-method}
\alias{[,bed.matrix,numeric,numeric,missing-method}
\alias{\%*\%,bed.matrix,matrix-method}
\alias{\%*\%,matrix,bed.matrix-method}
\alias{\%*\%,bed.matrix,vector-method}
\alias{\%*\%,vector,bed.matrix-method}
\alias{as.matrix,bed.matrix-method}
\alias{matrix,data.frameOrNULL,data.frameOrNULL-method}
\alias{standardize}
\alias{standardize,bed.matrix-method}
\alias{standardize<-}
\alias{standardize<-,bed.matrix-method}
\alias{coerce,bed.matrix,matrix-method}
\alias{coerce,matrix,bed.matrix-method}
\alias{coerce,bed.matrix,vector-method}
\alias{coerce,vector,bed.matrix-method}
\alias{dim,bed.matrix-method}
\alias{head,bed.matrix-method}
\alias{mu}
\alias{mu,bed.matrix-method}
\alias{mu<-}
\alias{mu<-,bed.matrix-method}
\alias{p}
\alias{p,bed.matrix-method}
\alias{p<-}
\alias{p<-,bed.matrix-method}
\alias{show,bed.matrix-method}
\alias{sigma}
\alias{sigma,bed.matrix-method}
\alias{sigma<-}
\alias{sigma<-,bed.matrix-method}
\alias{rbind}
\alias{rbind,bed.matrix-method}
\alias{cbind}
\alias{cbind,bed.matrix-method}

\title{Class \code{"bed.matrix"}}

\description{ S4 class for SNP genotype matrices }

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bed.matrix", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{ped}:}{ \code{data.frame} containing information for each individual: \code{famid} = Family ID, 
          \code{id} = Individual ID, \code{father} = Father ID, \code{mother} = Mother ID, \code{sex} = Sex and \code{pheno} = Phenotype. 
          Can also contain individuals statistic, for example: \code{N0}, \code{N1} and \code{N2} = Number of genotypes equal to 0, 1 and 2 respectively, 
          \code{NAs} = Number of missing genotypes, \code{callrate} = Individual callrate.}
    \item{\code{snps}:}{ \code{data.frame} containing information for each SNP: \code{chr} = Chromosome, \code{id} = SNP ID,
          \code{dist} = Genetic Distance, \code{pos} = Physical position, \code{A1} = Reference Allele, \code{A2} = Alternative Allele.
           Can also contain SNPs statistic, for example: \code{N0}, \code{N1} and \code{N2} = Number of genotypes equal to 0, 1 and 2 repectively,
           \code{NAs} = Number of missing genotypes, \code{callrate} = SNP callrate, \code{maf} = Minor allele frequency), \code{hz} = heterozygosity}
    \item{\code{bed}:}{ \code{externalptr} (pointing to the genotype matrix). }
    \item{\code{p}:}{ \code{vector} or \code{NULL} for allelic frequencies (allèle \code{A2}).}
    \item{\code{mu}:}{ \code{vector} or \code{NULL} for genotype means (usually \code{mu = 2*p}).}
    \item{\code{sigma}:}{ \code{vector} or \code{NULL} for genotypic standard deviation}
    \item{\code{standardize_p}:}{ \code{logical}. If \code{TRUE}, the genotype matrix is standardized using means \code{2*p} 
               and genotypic standard deviation \code{sqrt(2*p*(1-p))}}
    \item{\code{standardize_mu_sigma}:}{ \code{logical}. If \code{TRUE}, the genotype matrix is standardize using means
                \code{mu} and genotypic standard deviation \code{sigma}.}
  }
  For more details please check the vignette.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "bed.matrix", i = "numeric" or "logical" or "missing",}\cr\code{j = "numeric" or "logical" or "missing", drop = "missing")}: 
             \cr Extract a sub-matrix (a new \code{bed.matrix}). }
    \item{\%*\%}{\code{signature(x = "bed.matrix", y = "matrix" or "vector")}: 
             \cr Right matrix multiplication of the genotype matrix (possibly centered and reduced) with a \code{matrix} or a \code{vector}. }
    \item{\%*\%}{\code{signature(x = "matrix" or "vector", y = "bed.matrix")}: 
             \cr Left matrix multiplication of the genotype matrix with a \code{matrix} or a \code{vector}. }
    \item{as.matrix}{\code{signature(x = "bed.matrix")}: 
             \cr Convert a \code{bed.matrix} into a \code{matrix}.
	     The resulting matrix can be huge, use this method only for a small bed.matrix! }
    \item{standardize}{\code{signature(x = "bed.matrix")}: 
             \cr Get the standardize parameter of \code{bed.matrix}. Can be "none", "p" or "mu_sigma". }
    \item{standardize<-}{\code{signature(x = "bed.matrix")}:
             \cr Set the standardize parameter of a \code{bed.matrix}. }
    \item{dim}{\code{signature(x = "bed.matrix")}: 
             \cr Get the number of individuals (rows) and the number of SNPs (columns). }
    \item{head}{\code{signature(x = "bed.matrix")}: 
             \cr Print the head of the genotype matrix of a \code{bed.matrix} object. }
    \item{mu}{\code{signature(x = "bed.matrix")}: 
             \cr Get the \code{mu} slot of a \code{bed.matrix}. }
    \item{mu<-}{\code{signature(x = "bed.matrix")}: 
             \cr Set the \code{mu} slot of a \code{bed.matrix}. }
    \item{p}{\code{signature(x = "bed.matrix")}: 
             \cr Get the \code{p} slot of a \code{bed.matrix}. }
    \item{p<-}{\code{signature(x = "bed.matrix")}: 
             \cr Set the \code{p} slot of a \code{bed.matrix}. }
    \item{show}{\code{signature(object = "bed.matrix")}: 
             \cr Displays basic information about a \code{bed.matrix}. }
    \item{sigma}{\code{signature(x = "bed.matrix")}: 
             \cr Get the \code{sigma} slot of a \code{bed.matrix}. }
    \item{sigma<-}{\code{signature(x = "bed.matrix")}: 
             \cr Set the \code{sigma} slot of a \code{bed.matrix}. }
    \item{cbind}{\code{signature(... = "bed.matrix")}: 
             \cr Combine a sequence of \code{bed.matrix} by columns. }
    \item{rbind}{\code{signature(... = "bed.matrix")}: 
             \cr Combine a sequence of \code{bed.matrix} by rows. }
	 }
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{read.bed.matrix}}, \code{\link{write.bed.matrix}},
          \code{\link{set.stats}}, \code{\link{select.snps}}, \code{\link{select.inds}}, \code{\link{GRM}} }
\examples{
showClass("bed.matrix")

# Conversion example
data(LCT)
x1 <- as(LCT.gen, "bed.matrix")
x1
head(x1@ped)
head(x1@snps)

# the function as.bed.matrix is an alternative
x2 <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
x2
head(x2@ped)
head(x2@snps)
}

\keyword{classes}
