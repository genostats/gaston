\name{SNP.rm.duplicates}
\alias{SNP.rm.duplicates}

\title{ Remove duplicated SNPs }

\description{ Remove duplicated SNPs, taking into account possible genotype mismatches }
\usage{ SNP.rm.duplicates(x, by = "chr:pos", na.keep=TRUE, incomp.rm=TRUE) }

\arguments{
  \item{x}{ A bed.matrix }
  \item{by}{ The criterium used to determine duplicates }
  \item{na.keep}{ If \code{TRUE}, duplicated genotypes which are missing for at 
                  least one SNP are set to \code{NA}. }
  \item{incomp.rm}{ If \code{TRUE}, duplicated SNPs with allele incompatibility are
                  removed. Else, they are left as they are. }
}

\details{
Positions of duplicated SNPs are determined using \code{\link{SNP.duplicated}}
using parameter \code{by} (we recommend to use \code{"chr:pos"}, the default).

Then the function considers the possibility of alleles swaps or reference strand flips.
In case of allele incompatibility, the SNPs can be removed or kept as they are (\code{incomp.rm}
argument). In the contrary, only one copy is conserved. If there are genotype incompatibilities
between the duplicates for some individuals, these genotypes are set to \code{NA}.

Moreover the function takes special care of SNP with possible alleles \code{"0"}. 
This case occurs for monomorphic SNPs, when data are read from a \code{.ped} file; for
example, a whole column of \code{A A}'s will result in a SNP with alleles \code{"A"} and
\code{"0"}. If there's a duplicate of the SNP with a few, says, \code{A C}'s in it,
it will have alleles \code{"A"} and \code{"C"}. In that case, \code{\link{SNP.duplicated}} 
with \code{by = "chr:pos:alleles"} will not consider these SNPs as duplicates. 
}

\value{A bed.matrix without duplicated SNPs. }

\seealso{ \code{\link{SNP.match, SNP.duplicated}} }

\examples{
## Load simulations of 10 individuals with 7 duplicated SNPs
data(dupli)
x <- as.bed.matrix(dupli.gen, fam = dupli.ped, bim = dupli.bim)

# SNP 1 : No duplicated
# SNPs 2 and 3 : Two duplicated SNPs without alleles swap or reference strand flip
# SNPs 4 and 5 : Two duplicated SNPs with alleles swap
# SNPs 6 and 7 : Two duplicated SNPs with reference strand flip
# SNPs 8 and 9 : Two duplicated SNPs with alleles swap and reference strand flip
# SNPs 10 and 11 : Two duplicated SNPs with incompatible alleles
# SNPs 12 and 13 : Two duplicated SNPs including one monomorphic SNP (one allele to \code{"0"})
# SNPs 14, 15 and 16 : Three duplicated SNPs
# SNPs 17, 18 and 19 : Three duplicated SNPs with incompatible alleles
x@snps

# Missing values are included in genotypes data
as.matrix(x)

## By default (na.keep = TRUE), missing values are kept
as.matrix(SNP.rm.duplicates(x))

## Keep non missing genotypes
as.matrix(SNP.rm.duplicates(x, na.keep = FALSE))

## Keep SNPs with incompatible alleles
as.matrix(SNP.rm.duplicates(x, incomp.rm = FALSE))
}
