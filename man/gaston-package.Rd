\name{gaston-package}
\alias{gaston-package}
\alias{gaston}
\docType{package}
\title{ gaston }

\description{Manipulation of genetic data (SNPs), computation of Genetic Relationship Matrix, Linkage Disequilibrium, etc. 
Efficient algorithms for Linear Mixed Model (AIREML, diagonalisation trick).}

\section{Introducing gaston}{

  Gaston offers functions for efficient manipulation of 
  large genotype (SNP) matrices, and state-of-the-art implementation of algorithms
  to fit Linear Mixed Models, that can be used to compute heritability 
  estimates or to perform association tests.

  Thanks to the packages \code{\link[Rcpp:Rcpp-package]{Rcpp}}, 
  \code{\link[RcppParallel:RcppParallel-package]{RcppParallel}},
  \code{\link[RcppEigen:RcppEigen-package]{RcppEigen}}, gaston 
  functions are mainly written in C++. 

  Many functions are multithreaded;
  the number of threads can be setted through \code{RcppParallel}
  function \code{\link[RcppParallel:setThreadOptions]{setThreadOptions}}.
  It is advised to try several values for the number of threads, as 
  using too many threads might be conterproductive due to an important
  overhead.

  Some functions have a \code{verbose} argument, which controls the
  function verbosity. To mute all functions at once you can use 
  \code{options(gaston.verbose = FALSE)}.
}

\section{Genotype matrices}{
  An S4 class for genotype matrices is defined, named \code{\link{bed.matrix}}.
  Each row corresponds to an individual, and each column to a SNP. They can
  be read from files using \code{\link{read.bed.matrix}} 
  and saved using \code{\link{write.bed.matrix}}.  The function \code{\link{read.vcf}} reads
  VCF files.

  In first approach, a bed.matrix behaves as a "read-only" matrix containing only 
  0, 1, 2 and NAs, unless the genotypes are standardized (use \code{\link{standardize<-}}).
  They are stored in a compact form, each genotype being coded on 2 bits (hence
  4 genotypes per byte).

  Bed.matrices can be converted to numerical matrices with \code{\link{as.matrix}},
  and multiplied with numeric vectors or matrices with \code{\%*\%} (this 
  feature can be used e.g. to simulate quantitative phenotypes, see a basic example in the example
  section of \code{\link{association.test}}).

  It is possible to subset bed.matrices just as base matrices, writing e.g.
  \code{x[1:100,]} to extract the first 100 individuals, or \code{x[1:100,1000:1999]}
  for extract the SNPs 1000 to 1999 for these 100 individuals. The use of logical
  vectors for subsetting is allowed too. The functions 
  \code{\link{select.inds}} and \code{\link{select.snps}} can also be used for 
  subsetting with a nice syntax.

  Some basic descriptive statistics can be added to a bed.matrix with \code{\link{set.stats}} (since
  \code{gaston 1.4}, this function is called by default by all functions that create a bed.matrix, unless
  \code{options(gaston.auto.set.stats = FALSE)} was set.
  Hardy-Weinberg Equilibrium can be tested for all SNPs with \code{\link{set.hwe}}. 
}

\section{Crossproducts of standardized matrices}{

  If \eqn{X} is a standardized \eqn{n\times q}{n x q} genotype matrix, a Genetic Relationship Matrix 
  (GRM) of the individuals can be computed as
  \deqn{ GRM = {1\over q-1} XX’ }{GRM = XX’/(q-1)}
  where \eqn{q} is the number of SNPs. 
  This computation is done by the function \code{\link{GRM}}.  The eigen decomposition of the GRM produces
  the Principal Components (PC) of the data set. If needed, the loadings 
  corresponding to the PCs can be retrieved using \code{\link{bed.loadings}}.

  Doing the above crossproduct in the reverse order produces a moment estimate of the Linkage Disequilibrium:
  \deqn{ LD = {1\over n-1} X’X }{LD = X’X/(n-1)}
  where \eqn{n} is the number of individuals. This computation is done by the function
  \code{\link{LD}} (usually, only parts of the whole LD matrix is computed). This method is 
  also used by \code{\link{LD.thin}} to extract a set of SNPs in low linkage disequilibrium
  (it is often recommended to perform this operation before computing the GRM).
}

\section{Linear Mixed Models}{

  \code{\link{lmm.aireml}} is a function for linear mixed models parameter estimation
  and BLUP computations. 

  The model considered is of the form
  \deqn{ Y = X\beta + \omega_1 + \ldots + \omega_k + \varepsilon }{ Y = X beta + omega_1 + ... + omega_k + epsilon }
  with \eqn{ \omega_i \sim N(0,\tau_i K_i) }{omega_i ~ N(0, tau_i K_i)} for \eqn{ i \in 1, \dots,k }{i = 1, ..., k} and
  \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.

  Note that very often in genetics a mixed model is written as
  \deqn{ Y = X\beta + Zu + \varepsilon }{ Y = X beta + Zu + epsilon }
  with \eqn{Z} a standardized genotype matrix, and \eqn{u\sim N(0, \tau I_q)}{u ~ N(0, tau I_q)}. In that case, 
  denoting \eqn{\omega = Zu}{omega = Zu}, \eqn{\omega \sim N(0, \tau ZZ')}{omega ~ N(0,tau ZZ')}
  and letting \eqn{K=ZZ'} we get a mixed model of the previous form.

  When \eqn{k=1} in the above general model (only one random term \eqn{\omega}{omega}), the likelihood
  can be computed very efficiently using the eigen decomposition of 
  \eqn{K = \mathrm{var}(\omega)}{K = var(omega)}. This "diagonalization trick" 
  is used in \code{\link{lmm.diago.likelihood}} and \code{\link{lmm.diago}}, to compute 
  the likelihood and for parameter estimation, respectively. 
 
  Two small functions complete this set of functions: \code{\link{lmm.simu}}, to 
  simulate data under a linear mixed model, and \code{\link{random.pm}}, to generate
  random positive matrices. Both are used in examples and can be useful for data simulation.
}


\author{
Hervé Perdry and Claire Dandine-Roulland

Maintainer: Hervé Perdry
}

\keyword{Genetics}
\keyword{SNP}
\keyword{association study}
\keyword{linear mixed models}
