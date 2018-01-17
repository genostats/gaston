\name{SNP.match}
\alias{SNP.match}

\title{ SNP matching }

\description{ Returns a vector of the positions of (first) SNP matching of
              its first argument in its second. }
\usage{ SNP.match(x, table) }

\arguments{
  \item{x}{ A bed.matrix or a data.frame. }
  \item{table}{ A bed.matrix or a data.frame.}
}

\details{When \code{x} is a bed.matrix, the data.frame \code{x@bed} will be used; the 
same holds for \code{table}. The columns that will be taken in consideration
are \code{id}, \code{chr}, \code{pos}, \code{A1}, and \code{A2}. Not all columns
are mandatory (see below). All columns present in \code{x}
must be present in \code{table}.

The matching criterium depends on the columns present in \code{x}.
There are 5 possible criteria : (i) matching by chromosome and position, (ii)
matching by chromosome, position, and alleles, (iii) matching by id, 
(iv) matching by id, chromosome and position, and (v) matching by id, 
chromosome, position and alleles.

For each SNP in \code{x}, the function looks for the position of the first
matching SNP in \code{table}. If alleles are included in the matching criterium
(ie if allele columns \code{A1} and \code{A2} are present in x), the function 
also checks for SNP matching with swapped alleles (a SNP A/C would match a 
SNP C/A), or with reference strand flipped (i.e. a SNP A/C would match a SNP T/G) 
or both (a SNP A/C would match a SNP G/T).

This function should prove useful for data set merging.
}

\value{A named list with one or three members, depending on whether alleles are included
in the matching criterium.
  \item{index}{An integer vector giving the position of first match in \code{table}, or \code{NA} if there is no match}
  \item{swap}{A logical vector indicating whether the match is with swapped alleles}
  \item{flip}{A logical vector indicating whether the match is with flipped strand}
}

\seealso{ \code{\link{SNP.duplicated}} }

