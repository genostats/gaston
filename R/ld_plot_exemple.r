n <- 25; A <- matrix( runif(n**2,-1,1), nrow=n); A <- (A + t(A))/2; diag(A) <- 1
rownames(A) <- colnames(A) <- paste("rs", round(runif(n,0,1e7)),sep="" )

