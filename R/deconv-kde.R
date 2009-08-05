dkde.weights <- function(x, H, Sigma, ridge=FALSE, ...)
{
  require(kernlab)

  n <- nrow(x)
  d <- ncol(x)
  Lambda <- Sigma + H
  Omega <- Sigma + 2*H
  Q <- matrix(0, ncol=n, nrow=n)
  b <- rep(0, n)
  for (j in 1:n)
  {  
    xxj <- x - matrix(x[j,], nrow=n, ncol=d, byrow=TRUE)
    Q[j,] <- dmvnorm.mixt(xxj, mus=rep(0,d), Sigmas=2*Lambda, props=1)
    b[j] <- sum(dmvnorm.mixt(xxj, mus=rep(0,d), Sigmas=Omega, props=1))
  }
  Q <- Q/n^2
  b <- b/n^2
  if (ridge) Q <- Q + abs(min(Q))*1.1
  
  val <- ipop(c=-matrix(b, ncol=1), H=Q, A=matrix(1, ncol=n, nrow=1), b=n, r=0, l=matrix(0,ncol=1, nrow=n), u=matrix(n/10, ncol=1, nrow=n), ...)  

  return(primal(val))  
}

