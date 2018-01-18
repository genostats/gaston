random.pm <- function(n, values) {
  if(missing(values)) values <- grm.values(n)
  if(length(values) != n) stop("values should be of length n")
  Q <- random.ortho(n)
  K <- Q %*% (values * t(Q))
  return(list(K = K, eigen = list(values = values, vectors = Q)))
}

grm.values <- function(n) {
  kappa <- 1.63e-5
  gamma <- 0.39045 + n/13580
  delta <- 1/(1.1074 + n/175000)
  a <- 1/(1+exp( (gamma - qnorm(ppoints(n)))/delta ))
  mu.a <- mean(a)
  sd.a <- sd(a)
  rev( 1 + sqrt(kappa*n)*(a - mu.a)/sd.a )
}

random.ortho <- function(n)
  return(.Call('gg_random_ortho', PACKAGE = "gaston", n))

