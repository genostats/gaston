\name{SNP.duplicated}
\alias{SNP.duplicated}

\title{ Duplicated SNPs }

\description{ Determines which SNPs are duplicates of previous SNPs
              and returns their indices. }
\usage{ SNP.duplicated(x, by = "chr:pos") }

\arguments{
  \item{x}{ A bed.matrix or a data.frame with columns named \code{chr} and \code{pos}. }
  \item{by}{ The criterium used to determined if SNP is duplicated. Currently only one value possible. }
}

\details{When \code{by = "chr:pos"}, two SNPs are consider as duplicated if they have
same \code{chr} and \code{pos} values. Future versions will include other criteria.
}

\value{An integer vector of indices of SNPs which are duplicates of previously seen SNPs. }

\seealso{ \code{\link{SNP.match}} }

