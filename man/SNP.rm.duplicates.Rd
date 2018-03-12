\name{SNP.rm.duplicates}
\alias{SNP.rm.duplicates}

\title{ Remove duplicated SNPs }

\description{ Remove duplicated SNPs, taking into account possible genotype mismatches }
\usage{ SNP.rm.duplicates(x, by = "chr:pos", na.keep = TRUE, incomp.rm = TRUE) }

\arguments{
  \item{x}{ A bed.matrix }
  \item{by}{ The criterium used to determine duplicates }
  \item{na.keep}{ If \code{TRUE}, duplicated genotypes which are missing for at 
                  least one SNP are set to \code{NA}. }
  \item{incomp.rm}{ If \code{TRUE}, duplicated SNPs with allele incompatibility are
                    removed.}
}

\details{
Positions of duplicated SNPs are determined using \code{\link{SNP.duplicated}}
using parameter \code{by} (we recommend to use \code{"chr:pos"}, the default).

Then the function considers the possibility of alleles swaps or reference strand flips.
In case of allele incompatibility, the SNPs can be removed or not (according to \code{incomp.rm}
parameter).

When alleles can be matched, only one of the two SNPs is conserved. If there are 
genotype incompatibilities between the duplicates for some individuals, these genotypes are set 
to \code{NA}. The parameter \code{na.keep} settles the case of genotypes missing in one
of the SNPs.

Moreover the function takes special care of SNP with possible alleles \code{"0"}. 
This case occurs for monomorphic SNPs, when data are read from a \code{.ped} file; for
example, a whole column of \code{A A}'s will result in a SNP with alleles \code{"A"} and
\code{"0"}. If there's a duplicate of the SNP with a few, says, \code{A C}'s in it,
it will have alleles \code{"A"} and \code{"C"}. In that case, \code{\link{SNP.duplicated}} 
with \code{by = "chr:pos:alleles"} will not consider these SNPs as duplicates. 
}

\value{A bed.matrix without duplicated SNPs. }

\seealso{ \code{\link{SNP.match}}, \code{\link{SNP.duplicated}}, \code{\link{dupli}} }

\examples{
# Use example data of 10 individuals with 7 duplicated SNPs
data(dupli)
x <- as.bed.matrix(dupli.gen, fam = dupli.ped, bim = dupli.bim)

# There are any duplicated positions:
dupli.bim

x1 <- SNP.rm.duplicates(x)
# By default (na.keep = TRUE), as soon as the genotype is missing
# in one of the SNPs it is set to missing 
# (here looking at duplicated SNPs 2a and 2b)
as.matrix(x[,2:3])
as.matrix(x1[,2])

# With na.keep = FALSE 
x2 <- SNP.rm.duplicates(x, na.keep = FALSE)
as.matrix(x2[,2])

# Let's examinate SNP 3.a and 3.b (swapped alleles)
as.matrix(x[,4:5])
as.matrix(x1[,3])
as.matrix(x2[,3])

# and so on... (see also ?dupli)
}
