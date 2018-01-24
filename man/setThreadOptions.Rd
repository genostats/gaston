\name{setThreadOptions}
\alias{setThreadOptions}

\title{ Number of threads }

\description{ Set number of threads used in multi-threaded computations }
\usage{ setThreadOptions(...) }

\arguments{
  \item{...}{ Arguments similar to those of \code{RcppParallel::setThreadOptions} }
}

\details{
  Has it has been increasingly difficult to have packages using \code{RcppParallel} accepted
  on CRAN (see e.g. discussion at \url{https://github.com/RcppCore/RcppParallel/issues/36}) 
  Gaston has now is own multi-threading mechanism. It is only slightly less
  efficient than \code{RcppParallel} (by a few milliseconds due to thread creation overhead). 
  If you wish to enable \code{RcppParallel} in Gaston, define \code{GASTON_use_RcppParallel}
  to \code{true} in \code{src/Parallel.h} and recompile it.

  This function will work whether you compile with \code{RcppParallel} or not. Typical use
  is \code{ setTheadOptions(4) } to use 4 threads.
}

\seealso{ \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions}} }

