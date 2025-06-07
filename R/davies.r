### P value with Davies method
davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc",lambdas=as.double(lambda),noncentral=as.double(delta),df=as.integer(h),r=as.integer(r),sigma=as.double(sigma),q=as.double(q),lim=as.integer(lim),
             acc=as.double(acc),trace=as.double(rep(0,7)),ifault=as.integer(0),res=as.double(0), PACKAGE = "gaston")
  
  if (out$ifault==1) warning('In Davies method : Requested accuracy could not be obtained.')
  if (out$ifault==2) warning('In Davies method : Round-off error possibly significant.')
  if (out$ifault==3) stop('In Davies method : Invalid parameters.')
  if (out$ifault==4) stop('In Davies method : Unable to locate integration parameters.')
  
  out$res <- 1 - out$res
  if (out$res<0) out$res <- 0 
  return(out$res)
  
}

