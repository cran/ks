
###############################################################################
# Exact MISE for normal mixtures
###############################################################################

## nu, gamma.r, gamma.r2 written by Jose Chacon 10/2008

nu <- function(r,A)
{ ###Using the recursive formula provided in Kan (2008)
  A <- solve(A)
  ei <- eigen(A)$values
  tr <- numeric(r)
  for(p in 1:r)
    tr[p] <- sum(ei^p)

  nu <- 1
  if (r>=1)
  {
    for(p in 1:r)
    {
      a<-sum(tr[1:p]*rev(nu))/(2*p)
      nu<-c(nu,a)
    }
  } 
  return(factorial(r)*2^r*nu[r+1])   
}   

## gamma functional for MISE 
gamma.r <- function(mu, Sigma, d, r, Sd2r)
{
  Sigmainv <- chol2inv(chol(Sigma))
  w <- vec(K.pow(Sigmainv %*% Sigmainv, r)) %*% Sd2r
  v <- rep(0,length=d^(2*r))
  for(j in 0:r)
    v <- v + ((-1)^j*OF(2*j)*choose(2*r, 2*j))*(K.pow(mu,2*r-2*j)%x%K.pow(vec(Sigma),j))
  gamr <- (-1)^r*dmvnorm(mu,mean=rep(0,d),sigma=Sigma)*sum(w %*% v)
  
  return(gamr)
}

## gamma functional for AMISE 
gamma.r2 <- function(mu, Sigma, d, r, Sd2r4, H)
{
  Sigmainv <- chol2inv(chol(Sigma))
  w <- vec(K.pow(Sigmainv %*% Sigmainv, r)) %x% vec(K.pow(Sigmainv %*% H %*% Sigmainv, 2)) %*% Sd2r4
  ##w <- K.pow(vec(Sigmainv %*% Sigmainv), r) %x% K.pow(vec(Sigmainv %*% H %*% Sigmainv), 2) %*% Sd2r4 
  v <- rep(0,length=d^(2*r+4))
  for(j in 0:(r+2))
    v <- v+((-1)^j*OF(2*j)*choose(2*r+4, 2*j))*(K.pow(mu,2*r-2*j+4)%x%K.pow(vec(Sigma),j))
  
  gamr<-(-1)^r*dmvnorm(mu,mean=rep(0,d),sigma=Sigma)*sum(w %*% v)

  return(gamr)
}


###############################################################################
# Omega matrices (for exact MISE for normal mixtures)
#
# Parameters 
# mus - means
# Sigmas - variances
# k - number of mixture components
# a - subscript of Omega matrix
# H - bandwidth matrix
#
# Returns 
# Omega matrix
###############################################################################



omega <- function(mus, Sigmas, k, a, H, d, r, Sd2r)
{
  ## the (i,j) element of Omega matrix is dmvnorm(0, mu_i - mu_j,
  ## a*H + Sigma_i + Sigma_j)
 
  if (k == 1)
    omega.mat <- gamma.r(mu=rep(0,d),Sigma=a*H + 2*Sigmas, d=d, r=r, Sd2r=Sd2r)  ##dmvnorm(x=mus, mean=mus, sigma=a*H + 2*Sigmas)
  else
  {   
    if (is.matrix(mus)) d <- ncol(mus)
    else d <- length(mus)
    omega.mat <- matrix(0, nr=k, nc=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1):(i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1):(j*d),]
        muj <- mus[j,]    
        omega.mat[i,j] <- gamma.r(mu=mui-muj, Sigma=a*H + Sigmai + Sigmaj, d=d, r=r, Sd2r=Sd2r) ## dmvnorm(x=mui, mean=muj, sigma=a*H + Sigmai + Sigmaj)
      }
    }
  }
  
  return(omega.mat)
}


omega.1d <- function(mus, sigmas, k, a, h, d=1, r, Sd2r)
{
  ## the (i,j) element of Omega matrix is dmvnorm(0, mu_i - mu_j,
  ## a*H + sigma_i + sigma_j)

  H <- h^2
  Sigmas <- sigmas^2
  
  if (k == 1)
    omega.mat <- gamma.r(mu=0, Sigma=as.matrix(a*H + 2*Sigmas), d=d, r=r, Sd2r=Sd2r)  ##dmvnorm(x=mus, mean=mus, sigma=a*H + 2*Sigmas)
  else
  {   
    omega.mat <- matrix(0, nr=k, nc=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[i]
      mui <- mus[i]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[j]
        muj <- mus[j]    
        omega.mat[i,j] <- gamma.r(mu=mui-muj, Sigma=as.matrix(a*H + Sigmai + Sigmaj), d=d, r=r, Sd2r=Sd2r) ## dmvnorm(x=mui, mean=muj, sigma=a*H + Sigmai + Sigmaj)
      }
    }
  }
  
  return(omega.mat)
}




##############################################################################
# Exact MISE for normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# Props - vector of proportions of each mixture component 
# H - bandwidth matrix
# samp - sample size
#
# Returns
# Exact MISE for normal mixtures
###############################################################################

mise.mixt <- function(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
{
  if (!(missing(h)))
    return(mise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order)) 
  
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus)
  k <- length(props)
  r <- deriv.order
  Sd2r <- Sdr(d,2*r)

  ## formula is found in Wand & Jones (1993) and Chacon, Duong & Wand (2008)
  if (k == 1) 
  {
    mise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + 
        (1-1/samp)*omega(mus, Sigmas, 1, 2, H, d, r, Sd2r) -
        2*omega(mus, Sigmas, 1, 1, H, d, r, Sd2r) +
          omega(mus, Sigmas, 1, 0, H, d, r, Sd2r)
  }
  else
  {
    mise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      props %*% ((1-1/samp)*omega(mus, Sigmas, k, 2, H, d, r, Sd2r) - 
                 2*omega(mus, Sigmas, k, 1, H, d, r, Sd2r) + 
                 omega(mus, Sigmas, k, 0, H, d, r, Sd2r)) %*% props
  }
  return(drop(mise)) 
}

mise.mixt.1d <- function(h, mus, sigmas, props, samp, deriv.order=0)
{  
  d <- 1
  k <- length(props)
  r <- deriv.order
  Sd2r <- Sdr(d,2*r)
  H <- as.matrix(h^2)
  
  ## formula is found in Wand & Jones (1993) and Chacon, Duong & Wand (2008)
  if (k == 1) 
  {
    mise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + 
        (1-1/samp)*omega.1d(mus, sigmas, 1, 2, h, d, r, Sd2r) -
        2*omega.1d(mus, sigmas, 1, 1, h, d, r, Sd2r) +
          omega.1d(mus, sigmas, 1, 0, h, d, r, Sd2r)
  }
  else
  {
    mise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      props %*% ((1-1/samp)*omega.1d(mus, sigmas, k, 2, h, d, r, Sd2r) - 
                 2*omega.1d(mus, sigmas, k, 1, h, d, r, Sd2r) + 
                 omega.1d(mus, sigmas, k, 0, h, d, r, Sd2r)) %*% props
  }
  return(drop(mise)) 
}



 
###############################################################################
# Exact AMISE for bivariate normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# props - mixing proportions 
# H - bandwidth matrix
# samp - sample size
#
# Returns   
# Exact AMISE for normal mixtures
###############################################################################

amise.mixt <- function(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
{
  if (!(missing(h)))
    return(amise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))

  r <- deriv.order
  if (is.vector(mus)) {d <- length(mus); mus <- t(matrix(mus))}
  else d <- ncol(mus)
  k <- length(props)
 
  Sd2r4 <- Sdr(d,2*r+4)
  ##w <- Sd2r4%*%(K.pow(vec(diag(d)),r)%x%K.pow(vec(H),2))
  ##w <- as.vector(w)

  if (k == 1)
    omega.mat <- gamma.r2(mu=rep(0,d),Sigma=2*Sigmas, d=d, r=r, Sd2r4=Sd2r4, H=H)
  else
  {   
    omega.mat <- matrix(0, nr=k, nc=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1):(i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1):(j*d),]
        muj <- mus[j,]    
        omega.mat[i,j] <- gamma.r2(mu=mui-muj, Sigma= Sigmai + Sigmaj, d=d, r=r, Sd2r4=Sd2r4, H=H)
      }
    }
  }

  if (k == 1) {
    amise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + omega.mat/4
  }
  else {
    amise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      (props %*% omega.mat %*% props)/4
  }
 
  return(drop(amise))
}

amise.mixt.1d <- function(h, mus, sigmas, props, samp, deriv.order=0)
{
  d <- 1
  r <- deriv.order
  k <- length(props)
  H <- as.matrix(h^2)
  Sd2r4 <- Sdr(d,2*r+4)

  if (k == 1)
    omega.mat <- gamma.r2(mu=rep(0,d),Sigma=2*sigmas^2, d=d, r=r, Sd2r4=Sd2r4, H=H)
  else
  {   
    omega.mat <- matrix(0, nr=k, nc=k)
    for (i in 1:k)
    {
      Sigmai <- sigmas[i]^2
      mui <- mus[i]
      for (j in 1:k)
      {
        Sigmaj <- sigmas[j]^2
        muj <- mus[j]    
        omega.mat[i,j] <- gamma.r2(mu=mui-muj, Sigma= Sigmai + Sigmaj, d=d, r=r, Sd2r4=Sd2r4, H=H)
      }
    }
  }

  if (k == 1) {
    amise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) + omega.mat/4
  }
  else {
    amise <- 2^(-r)*nu(r,H)/(samp * (4 * pi)^(d/2) * sqrt(det(H))) +
      (props %*% omega.mat %*% props)/4
  }
 
  return(drop(amise))
}


amise.mixt.old <- function(H, mus, Sigmas, props, samp)
{ 
  if (is.vector(mus)) {d <- length(mus); mus <- t(matrix(mus))}
  else d <- ncol(mus)
  k <- length(props)
  Xi <- matrix(0, nr=k, nc=k)

  for (i in 1:k)
  {
    Sigmai <- Sigmas[((i-1)*d+1) : (i*d),]
    mui <- mus[i,]

    for (j in 1:k)
    {        
       Sigmaj <- Sigmas[((j-1)*d+1) : (j*d),]
       muj <- mus[j,]
       Aij <- chol2inv(chol(Sigmai + Sigmaj))
       Bij <- Aij %*% (diag(d) - 2*(mui - muj) %*%  t(mui - muj) %*% Aij)
       Cij <- Aij %*% (diag(d) - (mui - muj) %*%  t(mui - muj) %*% Aij)
    
       Xi[i,j] <- dmvnorm.mixt(x=mui, mus=muj, Sigmas=Sigmai+Sigmaj, props=1) *
                  (2*tr(H %*% Aij %*% H %*% Bij) + tr(H %*% Cij)^2)
    }  
  }
   
  amise <- 1/(samp *(4*pi)^(d/2)*sqrt(det(H)))+ 1/4*props %*% Xi %*% props

  return(drop(amise))
}


###############################################################################
# Lambda matrices (for exact AMISE for normal mixtures)
#
# Parameters 
# mus - means
# Sigmas - variances
# k - number of mixture components
# r - derivative (r1, r2)
#
# Returns 
# Lambda matrix
###############################################################################

lambda <- function(mus, Sigmas, k, r)
{
  ## the (i,j) element of Lambda matrix is d^r/ dx^r  dmvnorm(0, mu_i - mu_j,
  ## a*H + Sigma_i + Sigma_j)
    
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus)
  
  if (k == 1) 
    lambda.mat <- dmvnorm.deriv.2d(r=r, x=rep(0, length(mus)), Sigma=2*Sigmas)
  else
  {   
    if (is.matrix(mus)) d <- ncol(mus)
    else d <- length(mus)
    lambda.mat <- matrix(0, nr=k, nc=k)
    for (i in 1:k)
    {
      Sigmai <- Sigmas[((i-1)*d+1) : (i*d),]
      mui <- mus[i,]
      for (j in 1:k)
      {
        Sigmaj <- Sigmas[((j-1)*d+1) : (j*d),]
        muj <- mus[j,]    
        lambda.mat[i,j] <- dmvnorm.deriv.2d(r=r, x=mui-muj,Sigma=Sigmai+Sigmaj)
      }
    }
  }
  
  return(lambda.mat)
}


amise.mixt.2d <- function(H, mus, Sigmas, props, samp)
{  
  d <- ncol(Sigmas)
  k <- length(props)
  
  ##h1 <- sqrt(H[1,1])
  ##h2 <- sqrt(H[2,2])
  ##h12 <- H[1,2]

  ## formula is found in Wand & Jones (1993)
  if (k == 1) 
  {
   amise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) +
       1/4 *(lambda(mus, Sigmas, k, r=c(4,0))*H[1,1]^2 +
           4*lambda(mus, Sigmas, k, r=c(3,1))*H[1,1]*H[1,2] +  
           2*lambda(mus, Sigmas, k, r=c(2,2))*(H[1,1]*H[2,2] + 2*H[1,2]^2) + 
           4*lambda(mus, Sigmas, k, r=c(1,3))*H[2,2]*H[1,2]+    
             lambda(mus, Sigmas, k, r=c(0,4))*H[2,2]^2) 
  }
  else
  {
    amise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) +
      1/4 * props %*% 
          (  lambda(mus, Sigmas, k, r=c(4,0))*H[1,1]^2 +
           4*lambda(mus, Sigmas, k, r=c(3,1))*H[1,1]*H[1,2] +  
           2*lambda(mus, Sigmas, k, r=c(2,2))*(H[1,1]*H[2,2] + 2*H[1,2]^2) + 
           4*lambda(mus, Sigmas, k, r=c(1,3))*H[2,2]*H[1,2]+    
             lambda(mus, Sigmas, k, r=c(0,4))*H[2,2]^2) %*% props
  }
  
  return(drop(amise)) 
}


###############################################################################
# Finds the bandwidth matrix that minimises the MISE for normal mixtures
#
# Parameters
# mus - means
# Sigmas - variances
# props - vector of proportions of each mixture component 
# Hstart - initial bandwidth matrix
# samp - sample size
# full - 1 minimise over full bandwidth matrices
#      - 0 minimise over diagonal bandwidth matrices
# 
# Returns
# H_MISE
###############################################################################

hmise.mixt <- function(mus, sigmas, props, samp, hstart, deriv.order=0)
{
  r <- deriv.order
  d <- 1
  
  if (missing(hstart))
  {
    x <- rnorm.mixt(n=1000, mus=mus, sigmas=sigmas)  
    hstart <- sqrt((4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  mise.mixt.temp <- function(h)
  {  
    return(mise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optimize(f=mise.mixt.temp, interval=c(0, 10*hstart))
  hmise <- result$minimum
  return(hmise)
}

Hmise.mixt <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
 
  ## use normal reference estimate as initial condition
  if (missing(Hstart))
  {
    x <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- matrix.sqrt((4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  
  Hstart <- vech(Hstart)

  # input vech(H) into mise.mixt.temp because optim can only optimise
  # over vectors and not matrices
  mise.mixt.temp <- function(vechH)
  {  
    H <- invvech(vechH) %*% invvech(vechH)
    ## using H <- invvech(vechH) %*% invvech(vechH) ensures that H
    ## is positive definite
    
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optim(Hstart, mise.mixt.temp, method="BFGS")
  Hmise <- invvech(result$par) %*% invvech(result$par) 
  
  return(Hmise)
}   

Hmise.mixt.diag <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{   
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 

  if (missing(Hstart))
  {
    x <- rmvnorm.mixt(10000, mus, Sigmas, props)
    Hstart <- (4/(samp*(d + 2)))^(2/(d + 4)) * var(x)
  }  
  Hstart <- diag(matrix.sqrt(Hstart))

  mise.mixt.temp <- function(diagH)
  {  
    H <- diag(diagH) %*% diag(diagH)
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optim(Hstart, mise.mixt.temp, method = "BFGS")
  Hmise <- diag(result$par) %*% diag(result$par) 
  
  return(Hmise)
}   



###############################################################################
## Finds bandwidth matrix that minimises the AMISE for normal mixtures - 2-dim
##
## Parameters
## mus - means
## Sigmas - variances
## props - vector of proportions of each mixture component 
## Hstart - initial bandwidth matrix
## samp - sample size
## 
## Returns
## Bandwidth matrix that minimises AMISE
###############################################################################

hamise.mixt <- function(mus, sigmas, props, samp, hstart, deriv.order=0)
{
  r <- deriv.order
  d <- 1
 
  if (missing(hstart))
  {
    x <- rnorm.mixt(n=1000, mus=mus, sigmas=sigmas)  
    hstart <- sqrt((4/(samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  }
  amise.mixt.temp <- function(h)
  {  
    return(amise.mixt.1d(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp, deriv.order=deriv.order))
  }

  result <- optimize(f=amise.mixt.temp, interval=c(0, 10*hstart))
  hamise <- result$minimum
  
  return(hamise)
}


Hamise.mixt <- function(mus, Sigmas, props, samp, Hstart, deriv.order=0)
{
  r <- deriv.order
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 

  ## use explicit formula for single normal
  if (length(props)==1)
  {
    Hamise <- (4/ (samp*(d+2*r+2)))^(2/(d+2*r+4)) * Sigmas
  }
  else
  {  
    ## use normal reference estimate as initial condition
    if (missing(Hstart)) 
    {
      x <- rmvnorm.mixt(10000, mus, Sigmas, props)
      Hstart <- matrix.sqrt((4/ (samp*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
    }
    
    ## input vech(H) into mise.mixt.temp because optim can only optimise
    ## over vectors and not matrices    
    Hstart <- vech(Hstart)
    amise.mixt.temp <- function(vechH)
    {
      H <- invvech(vechH) %*% invvech(vechH)
      ## ensures that H is positive definite
      
      return(amise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=deriv.order))
    }
    
    result <- optim(Hstart, amise.mixt.temp, method="BFGS")
    Hamise <- invvech(result$par) %*% invvech(result$par) 
  }
  
  return(Hamise)
}   
  

###############################################################################
# ISE for normal mixtures (fixed KDE)
# 
# Parameters
# x - data values
# H - bandwidth matrix
# mus - matrix of means (each row is a vector of means from each component
#       density)
# Sigmas - matrix of covariance matrices (every d rows is a covariance matrix 
#          from each component density) 
# props - mixing proportions
# lower - vector of lower end points of rectangle
# upper - vector of upper end points of rectangle
# gridsize - vector of number of grid points
# stepsize - vector of step sizes
# Returns
# ISE 
###############################################################################

ise.mixt <- function(x, H, mus, Sigmas, props, h, sigmas, deriv.order=0)
{
  if (!(missing(h)))
    return(ise.mixt.1d(x=x, h=h, mus=mus, sigmas=sigmas, props=props, deriv.order=deriv.order))
  
  ##if (is.list(x))
  ##  return (ise.mixt.pc(x, H, mus, Sigmas, props, lower, upper, gridsize,
  ##                      stepsize))
  if (is.vector(x)) x <- matrix(x,nr=1)
  if (is.vector(mus)) mus <- matrix(mus, nr=length(props))

  d <- ncol(x)
  n <- nrow(x)
  M <- length(props)
  ise1 <- 0
  ise2 <- 0
  ise3 <- 0

  ## formula is found in thesis  
  if (d==2)
    ise1 <- dmvnorm.2d.sum(x=x, Sigma=2*H, inc=1)
  else if (d==3)
    ise1 <- dmvnorm.3d.sum(x=x, Sigma=2*H, inc=1)
  else if (d==4)
    ise1 <- dmvnorm.4d.sum(x=x, Sigma=2*H, inc=1)
  else if (d==5)
    ise1 <- dmvnorm.5d.sum(x=x, Sigma=2*H, inc=1)
  else if (d==6)
    ise1 <- dmvnorm.6d.sum(x=x, Sigma=2*H, inc=1)
  
  for (j in 1:M)
  {
    Sigmaj <- Sigmas[((j-1)*d+1) : (j*d),]
    ise2 <- ise2 + sum(props[j]*dmvnorm(x=x, mean=mus[j,], sigma=H + Sigmaj))
    
    for (i in 1:M)
    {
      Sigmai <- Sigmas[((i-1)*d+1) : (i*d),]
      ise3 <- ise3 + sum(props[i] * props[j] *
                         dmvnorm(x=mus[i,], mean=mus[j,], sigma=Sigmai+Sigmaj))
    }
  }  

  return (ise1/n^2 - 2*ise2/n + ise3)
}


ise.mixt.1d <- function(x, h, mus, sigmas, props, deriv.order=0)
{  
  d <- 1
  n <- length(x)
  M <- length(props)
  ise1 <- 0
  ise2 <- 0
  ise3 <- 0
  
  ise1 <- dnorm.sum(x=x, sigma=sqrt(2)*h, inc=1)
  
  for (j in 1:M)
  {
    sigmaj <- sigmas[j]
    ise2 <- ise2 + sum(props[j]*dnorm(x=x, mean=mus[j], sd=sqrt(h^2 + sigmaj^2)))
    
    for (i in 1:M)
    {
      sigmai <- sigmas[i]
      ise3 <- ise3 + sum(props[i]*props[j]*dnorm(x=mus[i], mean=mus[j], sd=sqrt(sigmai^2+sigmaj^2)))
    }
  }  

  return (ise1/n^2 - 2*ise2/n + ise3)
}


