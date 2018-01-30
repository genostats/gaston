\name{setThreadOptions}
\alias{setThreadOptions}
\alias{getNumThreads}

\title{ Number of threads }

\description{ Set ot get number of threads used in multi-threaded computations }
\usage{ setThreadOptions(numThreads = "auto", stackSize = "auto") 
        getNumThreads()}

\arguments{
  \item{numThreads}{ Number of threads to use }
  \item{stackSize}{ Stack size for threads }
}

\details{
  As it has been increasingly difficult to have packages using \code{RcppParallel} accepted
  on CRAN (see e.g. discussion at \url{https://github.com/RcppCore/RcppParallel/issues/36}) 
  Gaston uses now its own multi-threading mechanism. It is only slightly less
  efficient than \code{RcppParallel} (by at most a few milliseconds due to thread creation 
  overhead). 
  If you wish to enable \code{RcppParallel} in Gaston, define \code{GASTON_use_RcppParallel}
  to \code{true} in \code{src/Parallel.h} and recompile it.

  \code{setThreadOptions} will work whether you compile with \code{RcppParallel} or not. Typical use
  is \code{ setTheadOptions(4) } to use 4 threads. If you did not compile with \code{RcppParallel},
  parameter \code{stackSize} will be ignored.

  \code{getNumThreads} is not guaranteed to work if you compiled with \code{RcppParallel}.
  However, if the number of threads was set with this function it should send accurate results back.
}

\seealso{ \code{\link[RcppParallel:setThreadOptions]{RcppParallel::setThreadOptions}} }

