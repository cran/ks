
###############################################################################
# Exact MISE for normal mixtures
###############################################################################


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

omega <- function(mus, Sigmas, k, a, H)
{
  # the (i,j) element of Omega matrix is dmvnorm(0, mu_i - mu_j,
  # a*H + Sigma_i + Sigma_j)
  
  if (k == 1)
    omega.mat <- dmvnorm(x=mus, mean=mus, sigma=a*H + 2*Sigmas)
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
        omega.mat[i,j] <- dmvnorm(x=mui, mean=muj, sigma=a*H + Sigmai + Sigmaj)
      }
    }
  }
  
  return(omega.mat)
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
  # the (i,j) element of Lambda matrix is d^r/ dx^r  dmvnorm(0, mu_i - mu_j,
  # a*H + Sigma_i + Sigma_j)
    
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

mise.mixt <- function(H, mus, Sigmas, props, samp)
{  
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus)
  k <- length(props)

  ## formula is found in Wand & Jones (1993)
  if (k == 1) 
  {
    mise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) + 
      (1-1/samp)*omega(mus, Sigmas, 1, 2, H) -
        2*omega(mus, Sigmas, 1, 1, H) +
        omega(mus, Sigmas, 1, 0, H)
  }
  else
  {
    mise <- 1/(samp * (4*pi)^(d/2) * sqrt(det(H))) +
      props %*% ((1-1/samp)*omega(mus, Sigmas, k, 2, H) - 
                 2*omega(mus, Sigmas, k, 1, H) + 
                 omega(mus, Sigmas, k, 0, H)) %*% props
  }
  return(drop(mise)) 
}




 
###############################################################################
# Exact AMISE for bivariate normal mixtures - 2-dim
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



amise.mixt <- function(H, mus, Sigmas, props, samp)
{
  if (is.vector(mus)) {d <- length(mus); mus <- t(matrix(mus))}
  else d <- ncol(mus)
  k <- length(props)

  A <- matrix(0, nr=k, nc=k)
  B <- matrix(0, nr=k, nc=k)
  C <- matrix(0, nr=k, nc=k)
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

amise.mixt.2d <- function(H, mus, Sigmas, props, samp)
{  
  d <- ncol(Sigmas)
  k <- length(props)
  
  h1 <- sqrt(H[1,1])
  h2 <- sqrt(H[2,2])
  h12 <- H[1,2]

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
       
Hmise.mixt <- function(mus, Sigmas, props, samp, Hstart)
{   
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
  seed <- 8326

  # use normal reference estimate as initial condition
  set.seed(seed)
  x <- rmvnorm.mixt(1000, mus, Sigmas, props)
  if (missing(Hstart))
    Hstart <- matrix.sqrt((4/(samp*(d + 2)))^(2/(d + 4)) * var(x))
    
  Hstart <- vech(Hstart)

  # input vech(H) into mise.mixt.temp because optim can only optimise
  # over vectors and not matrices
  mise.mixt.temp <- function(vechH)
  {  
    H <- invvech(vechH) %*% invvech(vechH)
    # using H <- invvech(vechH) %*% invvech(vechH) ensures that H
    # is positive definite
    
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp))
  }

  result <- optim(Hstart, mise.mixt.temp, method="BFGS")
  Hmise <- invvech(result$par) %*% invvech(result$par) 
  
  return(Hmise)
}   

Hmise.mixt.diag <- function(mus, Sigmas, props, samp, Hstart)
{   
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
  seed <- 8326

  set.seed(seed)
  x <- rmvnorm.mixt(1000, mus, Sigmas, props)
  if (missing(Hstart))
    Hstart <- (4/(samp*(d + 2)))^(2/(d + 4)) * var(x)
    
  Hstart <- diag(matrix.sqrt(Hstart))

  mise.mixt.temp <- function(diagH)
  {  
    H <- diag(diagH) %*% diag(diagH)
    return(mise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp))
  }

  result <- optim(Hstart, mise.mixt.temp, method = "BFGS")
  Hmise <- diag(result$par) %*% diag(result$par) 
  
  return(Hmise)
}   



###############################################################################
# Finds bandwidth matrix that minimises the AMISE for normal mixtures - 2-dim
#
# Parameters
# mus - means
# Sigmas - variances
# props - vector of proportions of each mixture component 
# Hstart - initial bandwidth matrix
# samp - sample size
# 
# Returns
# Bandwidth matrix that minimises AMISE
###############################################################################
       
Hamise.mixt <- function(mus, Sigmas, props, samp, Hstart)
{   
  if (is.vector(mus)) d <- length(mus)
  else d <- ncol(mus) 
  seed <- 8326
  
  # use normal reference estimate as initial condition
  if (missing(Hstart)) 
  {
    set.seed(seed)
    x <- rmvnorm.mixt(1000, mus, Sigmas, props)
    Hstart <- matrix.sqrt((4/ (samp*(d + 2)))^(2/(d + 4)) * var(x))
  }
  
  # input vech(H) into mise.mixt.temp because optim can only optimise
  # over vectors and not matrices    
  Hstart <- vech(Hstart)
  amise.mixt.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    # ensures that H is positive definite
      
    return(amise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp))
  }
    
  result <- optim(Hstart, amise.mixt.temp, method="BFGS")
  Hamise <- invvech(result$par) %*% invvech(result$par) 
      
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

ise.mixt <- function(x, H, mus, Sigmas, props)
{  
  if (is.list(x))
    return (ise.mixt.pc(x, H, mus, Sigmas, props, lower, upper, gridsize,
                        stepsize))
  if (is.vector(x)) x <- matrix(x,nr=1)
  if (is.vector(mus)) mus <- matrix(mus, nr=length(props))
  
  d <- ncol(x)
  n <- nrow(x)
  M <- length(props)
  ise1 <- 0
  ise2 <- 0
  ise3 <- 0

  # formula is found in thesis  
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



###############################################################################
## ISE for t mixtures (numerical computation)
## 
## Parameters
## x - data values
## H - bandwidth matrices
## mus - matrix of means 
## Sigmas - matrix of covariance matrices 
## props - mixing proportions
## dfs - degrees of freedom
## lower, upper - lower and upper limits for integration
##
## Returns
## ISE  
###############################################################################

iset.mixt <- function(x, H, mus, Sigmas, dfs, props, lower,
                      upper, gridsize, stepsize=NULL) 
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))   
    stop("Proportions don't sum to one\n")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")  

  d <- ncol(Sigmas)
  if (missing(gridsize))
     gridsize <- rep(250, d) 

  ## pre-clustered KDE
  if (is.list(x))
  {
    x1 <- x
    x <- x$x
    n <-  nrow(x)
    ind.lab <- sort(unique(x1$ind))
    Hs <- numeric(0)
    for (i in 1:n)
    {
      clust <- which(x1$ind[i]==ind.lab)
      H1 <- H[((clust-1)*d+1):(clust*d),]
      Hs <- rbind(Hs, H1)
    }
  }
  ## fixed KDE
  else
  {
    n <-  nrow(x)
    Hs <- numeric(0)
    for (i in 1:n)
      Hs <- rbind(Hs, H) 
  }
  
  if(!is.null(gridsize))
  {  
    xx <- seq(lower[1], upper[1], length=gridsize[1])
    yy <- seq(lower[2], upper[2], length=gridsize[2])
  }
  else if (!is.null(stepsize))
  {
    xx <- seq(lower[1], upper[1], by=stepsize[1])
    yy <- seq(lower[2], upper[2], by=stepsize[2]) 
  } 
    
  xxyy <- permute(list(xx, yy))
  fhat <- dmvnorm.mixt(x=xxyy, mus=x, Sigma=Hs, props=rep(1/n,n))
  mixt <- dmvt.mixt(x=xxyy, mu=mus, Sigma=Sigmas, props=props, dfs=dfs)
  stepsize <- c(xx[1]-xx[2], yy[1]-yy[2])
  
  ise <- sum((fhat-mixt)^2*stepsize[1]*stepsize[2])

  return(ise)
}

###############################################################################
# Exact MISE for derivatives of normal density
###############################################################################


Hnr <- function(x, Hstart, deriv=0, pre="scale")   
{
  if(!is.matrix(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)

  pre1 <- substr(pre,1,2) 
  Sigma <- var(x)
  
  if (pre1=="sc") 
  {
     x.star <- pre.scale(x)
     S12 <- diag(sqrt(diag(var(x))))
  } 
  if (pre1=="sp") 
  {  
     x.star <- pre.sphere(x)
     S12 <- matrix.sqrt(var(x))
  }
     
  Sinv12 <- chol2inv(chol(S12))
  Sigma.star <- var(x.star)
 
  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma.star
  else
    Hstart <- Sinv12 %*% Hstart %*% Sinv12 
 
  if (deriv==0)
     H <- Hnr.dens(x=x.star, Hstart=Hstart)
  else if (deriv==1)
     H <- Hnr.grad(x=x.star, Hstart=Hstart)
  else if (deriv==2)
     H <- Hnr.curv(x=x.star, Hstart=Hstart)

      
  return(S12 %*% H %*% S12)
}


Hnr.diag <- function(x, Hstart, deriv=0, pre="scale")   
{
  if(!is.matrix(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)

  pre1 <- substr(pre,1,2) 
  if (pre1=="sp")
    stop("Using pre-sphering won't give a diagonal bandwidth matrix")

  Sigma <- var(x)
  
  if (pre1=="sc") 
  {
     x.star <- pre.scale(x)
     S12 <- diag(sqrt(diag(var(x))))
  }
  else if (pre1=="sp")
  {
     x.star <- pre.sphere(x)
     S12 <- matrix.sqrt(var(x))
  }
     
  Sinv12 <- chol2inv(chol(S12))
  Sigma.star <- var(x.star)
 
  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma.star
  else
    Hstart <- Sinv12 %*% Hstart %*% Sinv12 
 
  if (deriv==0)
     H <- Hnr.dens.diag(x=x.star, Hstart=Hstart)
  else if (deriv==1)
     H <- Hnr.grad.diag(x=x.star, Hstart=Hstart)
  else if (deriv==2)
     H <- Hnr.curv.diag(x=x.star, Hstart=Hstart)
      
  return(S12 %*% H %*% S12)
}


### Exact mise for density

mise.nr <- function(Sigma, d, n, H)
{  
  m.var <- n^(-1)*(4*pi)^(-d/2)*det(H)^(-1/2)
  m.bias1 <- (1-n^(-1)) * det(2*H + 2*Sigma)^(-1/2) 
  m.bias2 <- det(H + 2*Sigma)^(-1/2) 
  m.bias3 <- det(2*Sigma)^(-1/2) 

  return(m.var + (2*pi)^(-d/2)*(m.bias1 - 2*m.bias2 + m.bias3))
}


Hnr.dens <- function(x, Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma

  Hstart <- matrix.sqrt(Hstart)

  mise.nr.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(mise.nr(Sigma=Sigma, d=d, n=n, H=H))
  }

  result <- optim(vech(Hstart), mise.nr.temp, method="Nelder-Mead") 
  H <- invvech(result$par) %*% invvech(result$par)
  
  return(H)
}


Hnr.dens.diag <- function(x, Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma

  Hstart <- matrix.sqrt(Hstart)

  mise.nr.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(mise.nr(Sigma=Sigma, d=d, n=n, H=H))
  }

  result <- optim(diag(Hstart), mise.nr.temp, method="Nelder-Mead") 
  H <- diag(result$par) %*% diag(result$par)
  
  return(H)
}


### Exact MISE for gradient

mise.nr.grad <- function(Sigma, d, n, H)
{  
  mg.var <- 1/2*n^(-1)*(4*pi)^(-d/2)*det(H)^(-1/2) * chol2inv(chol(H))

  H2S2inv <- chol2inv(chol(2*H + 2*Sigma)) 
  HS2inv <- chol2inv(chol(H + 2*Sigma)) 
  S2inv <- chol2inv(chol(2*Sigma))
  mg.bias1 <- (1-n^(-1)) * det(2*H + 2*Sigma)^(-1/2) * H2S2inv
  mg.bias2 <- det(H + 2*Sigma)^(-1/2) * HS2inv
  mg.bias3 <- det(2*Sigma)^(-1/2) * S2inv

  return(tr(mg.var + (2*pi)^(-d/2)*(mg.bias1 - 2*mg.bias2 + mg.bias3)))
}


Hnr.grad <- function(x, Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma

  Hstart <- matrix.sqrt(Hstart)

  mise.nr.grad.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(mise.nr.grad(Sigma=Sigma, d=d, n=n, H=H))
  }

  result <- optim(vech(Hstart), mise.nr.grad.temp, method="Nelder-Mead") #, control=list(trace=2))
  H <- invvech(result$par) %*% invvech(result$par)
  
  return(H)
}

Hnr.grad.diag <- function(x, Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 6)) * Sigma

  Hstart <- matrix.sqrt(Hstart)
 
  mise.nr.grad.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    return(mise.nr.grad(Sigma=Sigma, d=d, n=n, H=H))  
  } 

  result <- optim(diag(Hstart), mise.nr.grad.temp, method="Nelder-Mead")
  H <- diag(result$par) %*% diag(result$par)
 
  return(H)
}

### Exact MISE for curvature

mise.nr.curv <- function(Sigma, d, n, H, R, D)
{  
  Hinv <- chol2inv(chol(H))
  if (missing(D))
  {
     Dd <- dupl(order=d)$d
     Dinv <- chol2inv(chol((t(Dd) %*% Dd)))
     D <- Dd %*% Dinv %*% Dinv %*% t(Dd)
  } 
  
  if (missing(R))
     R <- Rvec.D2(d=d)
  
  ##mc.var <- n^(-1)*det(H)^(-1/2) * (Hinv %x% Hinv) %*% R %*% D
  mc.var1 <- 2* tr((Hinv %x% Hinv) %*% D)
  mc.var2 <- t(vech(Hinv)) %*% vech(Hinv)  

  H2S2inv <- chol2inv(chol(2*H + 2*Sigma)) 
  HS2inv <- chol2inv(chol(H + 2*Sigma)) 
  S2inv <- chol2inv(chol(2*Sigma)) 

  mc.bias1 <- (1-n^(-1)) * det(2*H + 2*Sigma)^(-1/2) * (H2S2inv %x% H2S2inv) %*% D
  mc.bias2 <- det(H + 2*Sigma)^(-1/2) * (HS2inv %x% HS2inv) %*% D
  mc.bias3 <- det(2*Sigma)^(-1/2) * (S2inv %x% S2inv) %*% D

  return(1/4*(4*pi)^(-d/2)*n^(-1)*det(H)^(-1/2) * (mc.var1 + mc.var2) + tr(3*(2*pi)^(-d/2)*(mc.bias1 - 2*mc.bias2 + mc.bias3)))
}

Hnr.curv <- function(x, pre="scale", Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  Dd <- dupl(order=d)$d
  Dinv <- chol2inv(chol((t(Dd) %*% Dd)))
  D <- Dd %*% Dinv %*% Dinv %*% t(Dd)
  #R <- Rvec.D2(d=d)
 
  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma

  Hstart <- matrix.sqrt(Hstart)

  mise.nr.curv.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(mise.nr.curv(Sigma=Sigma, d=d, n=n, H=H, D=D))
  }

  result <- optim(vech(Hstart), mise.nr.curv.temp, method="Nelder-Mead")
  H <- invvech(result$par) %*% invvech(result$par)
  
  return(H)
}

Hnr.curv.diag <- function(x, Hstart)
{
  Sigma <- var(x)
  d <- ncol(x)
  n <- nrow(x)

  Dd <- dupl(order=d)$d
  Dinv <- chol2inv(chol((t(Dd) %*% Dd)))
  D <- Dd %*% Dinv %*% Dinv %*% t(Dd)
  #R <- Rvec.D2(d=d)

  if (missing(Hstart))
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * Sigma

  Hstart <- matrix.sqrt(Hstart)
   
  mise.nr.curv.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    return(mise.nr.curv(Sigma=Sigma, d=d, n=n, H=H, D=D))  
  } 

  result <- optim(diag(Hstart), mise.nr.curv.temp, method="Nelder-Mead")# control=list(trace=2))
  H <- diag(result$par) %*% diag(result$par)
 
  return(H)
}

### list element [[i]][[j]] is vector of indices for (i,j)-th element of x x^T
### e.g. x = (x1, x2)
### [[1]][[1]] = (2,0) <-> x1^2
### [[1]][[2]] = (1,1) <-> x1 * x2
### [[2]][[1]] = (1,1) <-> x1 * x2
### [[2]][[2]] = (0,2) <-> x2^2

xxt.index <- function(d)
{
  ind <- vector("list", d)
  for (i in 1:d)
     ind[[i]] <-  vector("list", d)
  
  for (i in 1:d)
   for (j in 1:d)
      ind[[i]][[j]] <- elem(i, d) + elem(j,d)
        
  return(ind)
}

### list element [[i]][[j]] is vector of indices for (i,j)-th element of
### x x^T kronecker x x^T 
### e.g. x = (x1, x2)
### [[1]][[1]] = (4,0) <-> x1^4
### [[1]][[2]] = (3,1) <-> x1^3 * x2
### [[1]][[3]] = (3,1) <-> x1^3 * x2
### [[1]][[4]] = (2,2) <-> x1^2 * x2^2 etc.

xxtxxt.index <- function(d)
{
  ind <- vector("list", d^2)
  for (i in 1:d^2)
    ind[[i]] <-  vector("list", d^2)
  
  xxt.ind <- xxt.index(d=d)
  
  for (i in 1:d)
    for (j in 1:d)
      for (k in 1:d)
        for (ell in 1:d)
          ind[[(i-1)*d + k]][[(j-1)*d + ell]] <- xxt.ind[[i]][[j]] + xxt.ind[[k]][[ell]]
  return(ind)   
}


### R (vec del^(2) f) where f is standard normal - used in mise.curv
### This is a d^2 x d^2 matrix

Rvec.D2 <- function(d)
{
  ind.coeff <- xxtxxt.index(d=d)
  mat.coeff <- matrix(0, ncol=d^2, nrow=d^2)
  
  for (i in 1:d^2)
    for (j in 1:d^2)
    {
      indij <- ind.coeff[[i]][[j]]
      ind.coeff[[i]][[j]][indij==4]  <- 3/4*(4*pi)^(-1/2)
      ind.coeff[[i]][[j]][indij==3]  <- 0
      ind.coeff[[i]][[j]][indij==2]  <- 1/2*(4*pi)^(-1/2)
      ind.coeff[[i]][[j]][indij==1]  <- 0
      ind.coeff[[i]][[j]][indij==0]  <- 1*(4*pi)^(-1/2)
      
      mat.coeff[i,j] <- prod(ind.coeff[[i]][[j]])
    }    
  
  return(mat.coeff)     
}


