#############################################################################
## Estimation of psi_r for normal scales
#############################################################################

psins.1d <- function(r, sigma)
{
  if (r %% 2 ==0)
    psins <- (-1)^(r/2)*factorial(r)/((2*sigma)^(r+1)*factorial(r/2)*pi^(1/2))
  else
    psins <- 0
    
  return(psins)  
}

psins.eta <- function(r, A, Sigma)
{
  if (r %% 2 ==0)
  {
    rr <- r/2 
    d <- ncol(Sigma)
    Sigmainv <- chol2inv(chol(Sigma))
    psins.val <- (4*pi)^(-d/2)*4^(-rr)*det(Sigma)^(-1/2)*nu.noncent(mu=rep(0,d), A=Sigmainv %*% A %*%Sigmainv,Sigma=-2*Sigma, r=r/2) 
  }
  else
    psins.val <- 0
  return(psins.val)
}

psins <- function(r, Sigma, deriv.vec=FALSE, Sdr.mat)
{
  d <- ncol(Sigma)
  if (deriv.vec)
  {
    return(drop(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), deriv.order=r, Sigma=2*Sigma, Sdr.mat=Sdr.mat)))
  }
  else
  {  
    if (d==2) return(psins.2d(r=r, Sigma=Sigma))
    ##if (d==3) return(psins.3d(r=r, Sigma=Sigma))
    ##if (d==4) return(psins.4d(r=r, Sigma=Sigma))
    ##if (d==5) return(psins.5d(r=r, Sigma=Sigma))
    ##if (d==6) return(psins.6d(r=r, Sigma=Sigma))
  }
}


psins.2d <- function(r, Sigma)
{
  A <- chol2inv(chol(2*Sigma))
  a11 <- A[1,1]
  a22 <- A[2,2]
  a12 <- A[1,2]

  
  if ((r[1] == 0) & (r[2] == 0))   
    psi <- 1
  
  ### second order psi functionals - normal score
  else if ((r[1]==2) & (r[2]==0))
    psi <- -a11
  else if ((r[1]==1) & (r[2]==1))
    psi <- -a12
  else if ((r[1]==0) & (r[2]==2))
    psi <- -a22
  
  ### fourth order psi functionals - normal score
  else if ((r[1]==4) & (r[2]==0))
    psi <-3*a11^2
  else if ((r[1]==3) & (r[2]==1))
    psi <-3*a11*a12
  else if ((r[1]==2) & (r[2]==2))
    psi <-2*a12^2 + a11*a22
  else if ((r[1]==1) & (r[2]==3))
    psi <-3*a12*a22
  else if ((r[1]==0) & (r[2]==4))
    psi <-3*a22^2
  
  ### sixth order psi functionals - normal score 
  else if ((r[1]==6) & (r[2]==0))
    psi <- -15*a11^3
  else if ((r[1]==5) & (r[2]==1))
    psi <- -15*a11^2*a12
  else if ((r[1]==4) & (r[2]==2))
    psi <- -3*a11*(4*a12^2 + a11*a22)
  else if ((r[1]==3) & (r[2]==3))
    psi <- -6*a12^3 - 9*a11*a12*a22
  else if ((r[1]==2) & (r[2]==4))
    psi <- -3*a22*(4*a12^2 + a11*a22)
  else if ((r[1]==1) & (r[2]==5))
    psi <- -15*a12*a22^2
  else if ((r[1]==0) & (r[2]==6))
    psi <- -15*a22^3

  ### eighth  order psi functionals - normal score 
  else if ((r[1]==8) & (r[2]==0))
    psi <- 105*a11^4
  else if ((r[1]==7) & (r[2]==1))
    psi <- 105*a11^3*a12
  else if ((r[1]==6) & (r[2]==2))
    psi <- 15*a11^2*(6*a12^2 + a11*a22)
  else if ((r[1]==5) & (r[2]==3))
    psi <- 15*a11*a12*(4*a12^2 + 3*a11*a22)
  else if ((r[1]==4) & (r[2]==4))
    psi <- 24*a12^4 + 72*a11*a12^2*a22 + 9*a11^2*a22^2
  else if ((r[1]==3) & (r[2]==5))
    psi <- 15*a12*a22*(4*a12^2 + 3*a11*a22)
  else if ((r[1]==2) & (r[2]==6))
    psi <- 15*a22^2*(6*a12^2 + a11*a22)
  else if ((r[1]==1) & (r[2]==7))
    psi <- 105*a12*a22^3
  else if ((r[1]==0) & (r[2]==8))
    psi <- 105*a22^4
  
  return(1/sqrt(2*pi)^2 * sqrt(det(A))*psi) 
}
  
