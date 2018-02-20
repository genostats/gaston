\name{SNP.duplicated.rm}
\alias{SNP.duplicated.rm}

\title{ Romve duplicated SNPs }

\description{ Remove duplicated SNPs considering mismatch between them. }
\usage{ SNP.duplicated.rm(x, by = "chr:pos") }

\arguments{
  \item{x}{ A bed.matrix }
  \item{by}{ The criterium used to determined if SNP is duplicated. }
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

For one individual, incompatible duplicated genotypes or genotypes which are missing for at least one SNP
are updtated to NA.
}

\value{A bed.matrix without duplicated SNPs. }

\seealso{ \code{\link{SNP.match, SNP.duplicated}} }

\examples{
set.seed(1)

## simulations of 10 individuals with 7 duplicated SNPs
bim <- data.frame(chr=1, id=paste0('SNP', 1:15, sep='_'), dist=0,
                  pos=c(rep(1:6, each=2), rep(7,3)), A1=NA, A2=NA)				  
fam <- data.frame(famid=1:10, id=1:10, father=0, mother=0, sex=0, pheno=0)

# Genotypes matrix
G <- matrix(NA, ncol=15, nrow=10)
## two identical SNPs
bim[1,c('A1','A2')] <- c('A', 'C')
bim[2,c('A1','A2')] <- c('A', 'C')
G[,1] <- G[,2] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[1,1] <- NA

## reference alllele reverse
bim[3,c('A1','A2')] <- c('A', 'C')
bim[4,c('A1','A2')] <- c('C', 'A')
G[,3] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,4] <- 2-G[,3]
G[2,4] <- NA

## reverse strand
bim[5,c('A1','A2')] <- c('A', 'C')
bim[6,c('A1','A2')] <- c('T', 'G')
G[,5] <- G[,6] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[3,5] <- G[4,6] <- NA

## reverse reference allele and strand
bim[7,c('A1','A2')] <- c('A', 'C')
bim[8,c('A1','A2')] <- c('G', 'T')
G[,7] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,8] <- 2-G[,7]
G[5,7] <- G[6,8] <- NA

## incompatible alleles
bim[9,c('A1','A2')] <- c('A', 'C')
bim[10,c('A1','A2')] <- c('A', 'T')
G[,9] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,10] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[7,9] <- NA

## one monomorphic SNP
bim[11,c('A1','A2')] <- c('A', 'C')
bim[12,c('A1','A2')] <- c('A', '0')
G[,11:12] <- 0
G[8:10, 11] <- c(1,1,NA)
G[8,12] <- NA

## more than 2 replications
bim[13,c('A1','A2')] <- c('A', 'C')
bim[14,c('A1','A2')] <- c('T', 'G')
bim[15,c('A1','A2')] <- c('G', 'T')
G[,13] <- G[,14] <- sample(c(0:2, NA), 10, replace=T, prob=c(0.45,0.30,0.20,0.05))
G[,15] <- 2-G[,13]
G[9,13] <- NA
G[10,15] <- NA
 
x <- as.bed.matrix(G, bim=bim, fam=fam)
xx <- SNP.duplicated.rm(x)

as.matrix(x)
as.matrix(xx)
}