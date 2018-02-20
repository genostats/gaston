\name{SNP.duplicated.rm}
\alias{SNP.duplicated.rm}

\title{ Romve duplicated SNPs }

\description{ Remove duplicated SNPs considering mismatch between them. }
\usage{ SNP.duplicated.rm(x, by = "chr:pos", na.keep=TRUE) }

\arguments{
  \item{x}{ A bed.matrix }
  \item{by}{ The criterium used to determined if SNP is duplicated. }
  \item{na.keep}{ If \code{TRUE}, duplicated genotypes which are missing for at least one SNP are updtated to NA. }
}

\details{
The columns that will be taken in consideration
Are \code{id}, \code{chr}, \code{pos}, \code{A1}, and \code{A2}. Not all columns
are mandatory, depending on the value of \code{by}.

The possible values for \code{by} are \code{"chr:pos"}, \code{"chr:pos:alleles"}, \code{"id"}, 
\code{"id:chr:pos"} and \code{"id:chr:pos:alleles"}. 
The default is \code{by = "chr:pos"}, which means that two SNPs are considered as duplicated if they have
same \code{chr} and \code{pos} values. 

This function considers the possibility of alleles swaps or reference strand flips. In case of
incompatibility, all duplicated SNPs are removed.

For one individual, incompatible duplicated genotypes are updtated to NA.
}

\value{A bed.matrix without duplicated SNPs. }

\seealso{ \code{\link{SNP.match, SNP.duplicated}} }

\examples{
set.seed(1)

## simulations of 10 individuals with 7 duplicated SNPs
bim <- data.frame(chr=1, id=paste0('SNP', 1:16, sep='_'), dist=0,
                  pos=c(rep(1:6, each=2), rep(7,3), 8), A1=NA, A2=NA)				  
fam <- data.frame(famid=1:10, id=1:10, father=0, mother=0, sex=0, pheno=0)

# Genotypes matrix
G <- matrix(NA, ncol=16, nrow=10)

## No duplicated
bim[1,c('A1','A2')] <- c('A', 'C')
G[,1] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))

## two identical SNPs
bim[2,c('A1','A2')] <- c('A', 'C')
bim[3,c('A1','A2')] <- c('A', 'C')
G[,2] <- G[,3] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[1,2] <- NA

## reference alllele reverse
bim[4,c('A1','A2')] <- c('A', 'C')
bim[5,c('A1','A2')] <- c('C', 'A')
G[,4] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,5] <- 2-G[,3]
G[2,5] <- NA

## reverse strand
bim[6,c('A1','A2')] <- c('A', 'C')
bim[7,c('A1','A2')] <- c('T', 'G')
G[,6] <- G[,7] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[3,7] <- G[4,7] <- NA

## reverse reference allele and strand
bim[8,c('A1','A2')] <- c('A', 'C')
bim[9,c('A1','A2')] <- c('G', 'T')
G[,8] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,9] <- 2-G[,8]
G[5,8] <- G[6,9] <- NA

## incompatible alleles
bim[10,c('A1','A2')] <- c('A', 'C')
bim[11,c('A1','A2')] <- c('A', 'T')
G[,10] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,11] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[7,10] <- NA

## one monomorphic SNP
bim[12,c('A1','A2')] <- c('A', 'C')
bim[13,c('A1','A2')] <- c('A', '0')
G[,12:13] <- 0
G[8:10, 12] <- c(1,1,NA)
G[8,13] <- NA

## more than 2 replications
bim[15,c('A1','A2')] <- c('A', 'C')
bim[15,c('A1','A2')] <- c('T', 'G')
bim[16,c('A1','A2')] <- c('G', 'T')
G[,14] <- G[,15] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,16] <- 2-G[,14]
G[9,14] <- NA
G[10,16] <- NA

 
x <- as.bed.matrix(G, bim=bim, fam=fam)
as.matrix(SNP.duplicated.rm(x))
as.matrix(SNP.duplicated.rm(x, na.keep=F))
}