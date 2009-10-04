dkde.weights <- function(x, H, Sigma, ridge=0, ...)
{
  require(kernlab)

  n <- nrow(x)
  d <- ncol(x)
  Lambda <- Sigma + H
  Omega <- Sigma + 2*H
  Qmat <- matrix(0, ncol=n, nrow=n)
  b <- rep(0, n)
  for (j in 1:n)
  {  
    xxj <- x - matrix(x[j,], nrow=n, ncol=d, byrow=TRUE)
    Qmat[j,] <- dmvnorm.mixt(xxj, mus=rep(0,d), Sigmas=2*Lambda, props=1)
    b[j] <- sum(dmvnorm.mixt(xxj, mus=rep(0,d), Sigmas=Omega, props=1))
  }

  if (ridge!=0)
  {
    eigval <- eigen(Qmat, only.values=TRUE, symmetric=TRUE)
    Qmat <- Qmat + diag(rep(max(c(abs(min(eigval$values))*ridge,1)), ncol(Qmat)))##max(c(abs(min(eigval$values)), 1))*ridge
  }
    
  Qmat <- Qmat/n^2
  b <- b/n^2
  
  val <- ipop(c=-matrix(b, ncol=1), H=Qmat, A=matrix(1, ncol=n, nrow=1), b=n, r=0, l=matrix(0,ncol=1, nrow=n), u=matrix(n, ncol=1, nrow=n), ...)  

  return(primal(val))  
}

