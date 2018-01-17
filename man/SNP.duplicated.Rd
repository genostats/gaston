\name{SNP.duplicated}
\alias{SNP.duplicated}

\title{ Duplicated SNPs }

\description{ Determines which SNPs are duplicates of previous SNPs
              and returns their indices. }
\usage{ SNP.duplicated(x, by = "chr:pos") }

\arguments{
  \item{x}{ A bed.matrix or a data.frame }
  \item{by}{ The criterium used to determined if SNP is duplicated. }
}

\details{
When \code{x} is a bed.matrix, the data.frame \code{x@bed} will be used.
The columns that will be taken in consideration
Are \code{id}, \code{chr}, \code{pos}, \code{A1}, and \code{A2}. Not all columns
are mandatory, depending on the value of \code{by}.

The possible values for \code{by} are \code{"chr:pos"}, \code{"chr:pos:alleles"}, \code{"id"}, 
\code{"id:chr:pos"} and \code{"id:chr:pos:alleles"}. 
The default is \code{by = "chr:pos"}, which means that two SNPs are considered as duplicated if they have
same \code{chr} and \code{pos} values. 

Currently, when using a criterium involving alleles, this function does not consider the possibility
of alleles swaps or reference strand flips.
}

\value{An integer vector of indices of SNPs which are duplicates of previously seen SNPs. }

\seealso{ \code{\link{SNP.match}} }

