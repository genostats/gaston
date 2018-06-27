\name{set.dist}
\alias{set.dist}

\title{ Set Genetic Distance }

\description{
  Returns an updated \code{\link{bed.matrix}} with genetic 
  distances in centimorgan computed from the variant positions
}
\usage{ set.dist(x, map, verbose = getOption("gaston.verbose", TRUE)) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{map}{ The genetic map, given by a list of data frames (see Details) }
 \item{verbose}{ If \code{TRUE}, display information on the function actions }
}

\details{
  A map is a list of data frames, with names corresponding to chromosomes.
  Each of these data frames must have columns \code{pos} and \code{dist} corresponding
  to positions in bp and cM, respectively.

  Such maps are too large to be included in a CRAN package. You can get two genetic
  maps for the Human Genome (build 36 and 37) in the package \code{HumanGeneticMap}  
  on GitHub.

  To install this package, run

\code{ install.packages("HumanGeneticMap", repos="https://genostats.github.io/R/") }

  You can then use this function with \code{set.dist(x, HumanGeneticMap::genetic.map.b36)}
  for example, for positions on the build 36. Use \code{map = HumanGeneticMap::genetic.map.b37}) 
  for the build 37.

}

\value{
  A \code{\link{bed.matrix}} similar to \code{x}, with updated values in \code{x@snps$dist}.
}

