

###############################################################################
## Estimate g_AMSE pilot bandwidths for even orders - 2-dim
##
## Parameters
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
##
## Returns
## g_AMSE pilot bandwidths for even orders
###############################################################################

gamse.even.2d <- function(r, n, psi1, psi2)
{
  d <- 2
  num <- -2 * dmvnorm.deriv(x=c(0,0), deriv.order=r, Sigma=diag(2))
  den <- (psi1 + psi2) * n   
  g.amse <- (num/den)^(1/(2 + d + sum(r)))
  
  return(g.amse)
}

###############################################################################
## Estimate g_AMSE pilot bandwidths for odd orders - 2-dim
##
## Parameters
##
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
## psi00 - psi_(0,0)
## RK - R(K^(r))
## 
## Returns
## g_AMSE pilot bandwidths for odd orders
###############################################################################

gamse.odd.2d <- function(r, n, psi1, psi2, psi00, RK)
{  
    d <- 2
    num <- 2 * psi00 * (2 * sum(r) + d) * RK
    den <- (psi1 + psi2)^2 * n^2
    g.amse <- (num/den)^(1/(2*sum(r) + d + 4))

    return(g.amse)
}


###############################################################################
## Estimate g_SAMSE pilot bandwidth - 2- to 6-dim 
##
## Parameters
## Sigma.star - scaled variance matrix
## n - sample size
##
## Returns
## g_SAMSE pilot bandwidth
###############################################################################

gsamse.nd <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- ncol(Sigma.star)
  K <- numeric(); psi <- numeric()

  derivt4 <- deriv.list(d=d, r=4)
  derivt6 <- deriv.list(d=d, r=6)
 
  ## 4th order g_SAMSE
  if (modr == 4)
  {
    for (i in 1:nrow(derivt4))
    {
      r <- derivt4[i,]
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv(x=rep(0,d), deriv.order=r, Sigma=diag(d)))
        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
            A3psi <- A3psi + psins(r=r+2*elem(j,d), Sigma=Sigma.star)
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
        }
        psi <- c(psi, A3psi)    
      }
    }
  }
  ## 6th order g_SAMSE
  else if (modr==6)
  {
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv(x=rep(0,d), deriv.order=r, Sigma=diag(d)))
       
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins(r=r+2*elem(j,d), Sigma=Sigma.star)
        psi <- c(psi, A3psi)
      }
    }
  }
  
  ## see thesis for formula
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}

##############################################################################
## Estimate psi functionals for bivariate data using 1-stage plug-in - 2-dim
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
##
## Returns
## estimated psi functionals
###############################################################################

psifun1.2d <- function(x.star, pilot="samse", binned, bin.par)
{
  d <- 2
  derivt4 <- deriv.list(d=d, r=4) ##cbind((2*d) - 0:(2*d), 0:(2*d)) 
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  RK31 <- 15/(64*pi)
  psi00 <- psins(r=c(0,0), Sigma=S.star) 
  psihat.star <- vector()
  g.star <- vector()
  
  ## pilots are based on 4th order derivatives
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- rep(gsamse.nd(S.star, n, 4), nrow(derivt4))
 
  ## compute 5 different pilots for AMSE
  else if (pilot=="amse")
    for (k in 1:nrow(derivt4))
    { 
      r <- derivt4[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)

      ## odd order
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK31)
      
      ## even order
      else
        g.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }
  
  if (!binned)
    x.star.diff <- differences(x.star, upper=FALSE)
   
  for (k in 1:nrow(derivt4))
  {
    r <- derivt4[k,]
    G.star <- g.star[k]^2 * diag(2)
    if (binned)
      psihat.star[k] <- kfe(bin.par=bin.par, G=G.star, r=r, binned=TRUE)
    else 
      psihat.star[k] <- kfe.scalar(x=x.star.diff, r=r, g=g.star[k], diff=TRUE) 
  }
 
  return(psihat.star)
}


###############################################################################
# Estimate psi functionals for bivariate data using 2-stage plug-in - 2-dim
#
# Parameters
# x - pre-transformed data points
# pilot - "amse" - different AMSE pilot
#       - "samse" - SAMSE pilot
# Returns
# estimated psi functionals
###############################################################################

psifun2.2d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 2
  derivt4 <- deriv.list(d=d, r=4)
  derivt6 <- deriv.list(d=d, r=6)
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  RK31 <- 15/(64*pi)
  RK51 <- 945/(256*pi)
  RK33 <- 225/(256*pi)
  psi00 <- psins(r=c(0,0), Sigma=S.star) 

  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- rep(gsamse.nd(S.star, n, 6), nrow(derivt6))  
    
  ## compute different pilots for AMSE
  else if (pilot=="amse")
  {       
    for (k in 1:nrow(derivt6))
    {
      r <- derivt6[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)
      if (prod(r) == 5)
        g6.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK51)
      else if (prod(r) == 9)
        g6.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK33) 
      else  
        g6.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }
  }

  if (!binned) x.star.diff <- differences(x.star, upper=FALSE)
  
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star[k]^2 * diag(d)

    if (binned)
      psihat6.star[k] <- kfe(bin.par=bin.par, G=G6.star, r=r, binned=TRUE)
    else
      psihat6.star[k] <- kfe.scalar(x=x.star.diff, r=r, g=g6.star[k], diff=TRUE)
  }

 
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- rep(gsamse.nd(S.star, n, 4, nstage=2, psihat=psihat6.star), nrow(derivt4))  
  else if (pilot=="amse")
    for (k in 1:nrow(derivt4))
    {
      r <- derivt4[k,]
      psi1 <- psihat6.star[7 - (r + 2*elem(1,2))[1]]
      psi2 <- psihat6.star[7 - (r + 2*elem(2,2))[1]]
      
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK31)
      else
        g.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }
  
  for (k in 1:nrow(derivt4))
  {
    r <- derivt4[k,]
    G.star <- g.star[k]^2 * diag(2)

    if (binned)
      psihat.star[k] <- kfe(bin.par=bin.par, G=G.star, r=r, binned=TRUE)
    else 
      psihat.star[k] <- kfe.scalar(x=x.star.diff, r=r, g=g.star[k], diff=TRUE)
  }

  return(psihat.star)
}


###############################################################################
## Estimate psi functionals for 3-variate data using 1-stage plug-in - 3-dim
##
## Parameters
## x.star - pre-transformed data points
## pilot - "samse" = optimal SAMSE pilot bandwidth
## Returns
## estimated psi functionals
###############################################################################

psifun1.nd <- function(x.star, d, pilot="samse", binned, bin.par)
{ 
  derivt <- Psi4.list(d)$psi
  derivt4 <- deriv.list(d, r=4)
  S.star <- var(x.star)
  n <- nrow(x.star)

  psihat.star <- rep(0, length=nrow(derivt))
  g.star <- vector()

  if (!binned) x.star.diff <- differences(x.star, upper=FALSE)
   
  ## compute 1 pilot for SAMSE
  g.star <- gsamse.nd(S.star, n, 4, nstage=1)
  G.star <- g.star^2 * diag(d)
 
  for (k in 1:nrow(derivt4))
  {
    r <- derivt4[k,]
    kind <- which.mat(r, derivt)
    if (binned)
      psihat.star[kind] <- kfe(bin.par=bin.par, G=G.star, r=r, binned=TRUE)
    else 
      psihat.star[kind] <- kfe.scalar(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }

  return(psihat.star)
}



###############################################################################
# Estimate psi functionals for 3-variate data using 2-stage plug-in - 3-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
#
# Returns
# estimated psi functionals
###############################################################################

psifun2.nd <- function(x.star, d, pilot="samse", binned, bin.par)
{
  derivt <- Psi4.list(d)$psi
  derivt4 <- deriv.list(d, r=4)
  derivt6 <- deriv.list(d, r=6)
 
  S.star <- var(x.star)
  n <- nrow(x.star)

 
  if (!binned)
    x.star.diff <- differences(x.star, upper=FALSE)

  psihat6.star <- vector()
  g6.star <- vector()
  psihat.list.star <- vector()
  psihat.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
   
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.nd(Sigma.star=S.star, n=n, modr=6)
  G6.star <- g6.star^2 * diag(d) 
  
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    if (binned)
      psihat6.star[k] <- kfe(bin.par=bin.par, G=G6.star, r=r, binned=TRUE)
    else 
      psihat6.star[k] <- kfe.scalar(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.nd(S.star, n, 4, nstage=2, psihat=psihat6.star) 

  G.star <- g.star^2 * diag(d)
  
  for (k in 1:nrow(derivt4))
  {
    r <- derivt4[k,]
    kind <- which.mat(r, derivt)
    if (binned)
      psihat.star[kind] <- kfe(bin.par=bin.par, G=G.star, r=r, binned=TRUE)
    else 
      psihat.star[kind] <- kfe.scalar(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }  

  return(psihat.star)
}


#############################################################################
## Estimate psi functionals for 6-variate data using 1-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
#############################################################################

psifun1.unconstr.nd <- function(x, Sd4, Sd6, rel.tol=10^-10)
{
  n <- nrow(x)
  d <- ncol(x)
  S <- var(x)
  
  nlim <- 1e4
  upper <- TRUE 

  ## stage 1 of plug-in
  G4 <-(2^(d/2+3)/((d+4)*n))^(2/(d+8))*S
  vecPsi4 <- vecPsir(x=x, Sdr=Sd4, Gr=G4, r=4, upper=upper, nlim=nlim)
  
  return (vecPsi4)
}


#############################################################################
## Estimate psi functionals for 6-variate data using 2-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
############################################################################

psifun2.unconstr.nd <- function(x, Sd4, Sd6, rel.tol=10^-10)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  Hstart <- (4/(d+2))^(2/(d+4))*n^(-2/(d+4))*S
  Hstart <- matrix.sqrt(Hstart)
  nlim <- 1e4
 
  ## matrix of pairwise differences
  upper <- TRUE
  difs <- differences(x, upper=upper)
 
  ## constants for normal reference
  D4phi0 <- D4L0(d=d, Sd4=Sd4)
  Id1 <- diag(d)
  vId <- vec(Id1)

  ## stage 1 of plug-in
  G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
  G612 <- matrix.sqrt(G6)
 
  vecPsi6 <- vecPsir(x=x, Sdr=Sd6, Gr=G6, r=6, upper=upper, nlim=nlim)
   
  Id4 <- diag(d^4)
  Id2 <- diag(d^2)
  Kdd2 <- K.mat(m=d,n=d^2)
  Psi6 <- (Id1%x%Kdd2%x%Id2)%*%vecPsi6
  
  ## asymptotic squared bias for r = 4
  AB2r4<-function(vechG){
    r <- 4
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=r)%*%D4phi0)+
      (1/2)*(t(vec(G))%x%Id4)%*%Psi6
    return (sum(AB^2))
  }
  
  res <- optim(vech(Hstart),AB2r4, control=list(reltol=rel.tol))
  V4 <- res$value
  G4 <- res$par
  G4 <- invvech(G4)%*%invvech(G4)
 
  ## stage 2 of plug-in

  vecPsi4 <- vecPsir(x=x, Sdr=Sd4, Gr=G4, r=4, upper=upper, nlim=nlim)
  
  return (vecPsi4)
}


############################################################################
## Psi_4 matrix of 4th order psi functionals used in AMISE - 2 to 6 dim
##
## Parameters
## x - data points
## nstage - number of plug-in stages (1 or 2)
## pilot - "amse" - different AMSE pilot
##       - "samse" - SAMSE pilot
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## matrix of psi functionals
############################################################################

psimat.2d <- function(x.star, nstage=1, pilot="samse", binned, bin.par)
{
  if (nstage==1)
    psi.fun <- psifun1.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  else if (nstage==2)
    psi.fun <- psifun2.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)

  psi40 <- psi.fun[1] 
  psi31 <- psi.fun[2] 
  psi22 <- psi.fun[3] 
  psi13 <- psi.fun[4] 
  psi04 <- psi.fun[5] 

  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13, psi04)

  return(matrix(coeff * psi.fun, nc=3, nr=3))
}

psimat.nd <- function(x.star, d, nstage=1, pilot="samse", binned, bin.par)
{
  if (nstage==1)
    psi.fun <- psifun1.nd(x.star, d=d, pilot=pilot, binned=binned, bin.par=bin.par)
  else if (nstage==2)
    psi.fun <- psifun2.nd(x.star, d=d, pilot=pilot, binned=binned, bin.par=bin.par)

  coeff <- Psi4.list(d)$coeff
  
  return(matrix(coeff * psi.fun, nc=d*(d+1)/2, nr=d*(d+1)/2))
}

###  unconstrained pilot selectors

psimat.unconstr.nd <- function(x, nstage=1, Sd4, Sd6)
{
  if (nstage==1)
    psi.fun <- psifun1.unconstr.nd(x=x, Sd4=Sd4, Sd6=Sd6)
  else if (nstage==2)
    psi.fun <- psifun2.unconstr.nd(x=x, Sd4=Sd4, Sd6=Sd6)

  return(invvec(psi.fun))
}
  

#############################################################################
# Plug-in bandwidth selectors
#############################################################################

    
############################################################################
## Computes plug-in full bandwidth matrix - 2 to 6 dim
##
## Parameters
## x - data points
## Hstart - initial value for minimisation
## nstage - number of plug-in stages (1 or 2)
## pilot - "amse" - different AMSE pilot
##       - "samse" - SAMSE pilot
##       - "unconstr" - unconstrained pilot
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## Plug-in full bandwidth matrix
###############################################################################

hpi <- function(x, nstage=2, binned=TRUE, bgridsize)
{
  ## 1-d selector is taken from KernSmooth's dpik
  
  if (missing(bgridsize)) bgridsize <- 401
  return(dpik(x=x, level=nstage, gridsize=bgridsize))
}

Hpi <- function(x, nstage=2, pilot="samse", pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hpi", kfold=kfold, random=FALSE, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE, amise=FALSE))
    else
      return(Hkfold(x=x, selector="Hpi", kfold=kfold, random=FALSE, Hstart=Hstart, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE, amise=FALSE))
  }

  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)

  if(!is.matrix(x)) x <- as.matrix(x)

  if (substr(pre,1,2)=="sc")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (substr(pre,1,2)=="sp")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }
  
  if (substr(pilot,1,1)=="a")
    pilot <- "amse"
  else if (substr(pilot,1,1)=="s")
    pilot <- "samse"
  else if (substr(pilot,1,1)=="u")
    pilot <- "unconstr"
           
  if (pilot=="amse" & d>2)
    stop("SAMSE pilot selectors are better for higher dimensions")

  if (pilot=="unconstr" & d>=6)
    stop("Uconstrained pilots not implemented yet for 6-dim data")
  
  if (missing(bgridsize) & binned) bgridsize <- default.gridsize(d)
  if (d > 4) binned <- FALSE
  
  if (binned)
  {
    ## linear binning
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    ##bin.par <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max)))

    ## for large samples, take subset for pilot estimation
    nsub <- n ## min(n, 1e4)
    x.star.sub <- x.star[sample(1:n, size=nsub),] 
    bin.par.sub <- binning(x=x.star.sub, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
  }
  else 
    x.star.sub <- x.star
 
  ## psi4.mat is on pre-transformed data scale
  if (pilot!="unconstr")
  {
    if (d==2)
      psi4.mat <- psimat.2d(x.star.sub, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par.sub)
    else 
      psi4.mat <- psimat.nd(x.star.sub, d=d, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par.sub)
  }
  else
  {
    ### symmetriser matrices for unconstrained pilot selectors
    Sd4 <- Sdr(d=d, r=4)
    Sd6 <- Sdr(d=d, r=6)
    psi4.mat <- psimat.unconstr.nd(x=x, nstage=nstage, Sd4=Sd4, Sd6=Sd6)
  }

  if (pilot=="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x)
    
    Hstart <- matrix.sqrt(Hstart)

    ## PI is estimate of AMISE
    pi.temp.unconstr <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vec(H)) %*% psi4.mat %*% vec(H)
      return(drop(pi.temp))
    }

    ## psi4.mat always a zero eigen-value since it has repeated rows 
    result <- optim(vech(Hstart), pi.temp.unconstr, method="BFGS")
    H <- invvech(result$par) %*% invvech(result$par)
  }
  else if (pilot!="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else
      Hstart <- Sinv12 %*% Hstart %*% Sinv12
    
    Hstart <- matrix.sqrt(Hstart)

    ## PI is estimate of AMISE
    pi.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vech(H)) %*% psi4.mat %*% vech(H)
      return(drop(pi.temp)) 
    }
    
    ## check that Psi_4 is positive definite  
    if (prod(eigen(psi4.mat)$val > 0) == 1)
    {
      result <- optim(vech(Hstart), pi.temp, method="BFGS")
      H <- invvech(result$par) %*% invvech(result$par)
    }
    else
    { 
      cat("Psi matrix not positive definite\n")
      H <- matrix(NA, nc=d, nr=d) 
    }    
    ## back-transform
    H <- S12 %*% H %*% S12
  }

  if (!amise)
    return(H)
  else
    return(list(H = H, PI=result$value))
 
}     



###############################################################################
# Computes plug-in diagonal bandwidth matrix for 2 to 6-dim
#
# Parameters
# x - data points
# nstage - number of plug-in stages (1 or 2)
# pre - "scale" - pre-scaled data
#     - "sphere"- pre-sphered data 
#
# Returns
# Plug-in diagonal bandwidth matrix
###############################################################################


Hpi.diag <- function(x, nstage=2, pilot="amse", pre="scale", Hstart, binned=FALSE, bgridsize, kfold=1)
{
  if(!is.matrix(x)) x <- as.matrix(x)

  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hpi.diag", kfold=kfold, random=FALSE, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE))
    else
      return(Hkfold(x=x, selector="Hpi.diag", kfold=kfold, random=FALSE, Hstart=Hstart, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE))
  }
  
  if (substr(pre,1,2)=="sc")
    x.star <- pre.scale(x)
  else if (substr(pre,1,2)=="sp")
    x.star <- pre.sphere(x)

  if (substr(pre,1,2)=="sp")
    stop("Using pre-sphering won't give diagonal bandwidth matrix\n")

  if (substr(pilot,1,1)=="a")
    pilot <- "amse"
  else if (substr(pilot,1,1)=="s")
    pilot <- "samse"
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  s1 <- sd(x[,1])
  s2 <- sd(x[,2])

  if (substr(pre,1,2)=="sc") S12 <- diag(sqrt(diag(var(x))))
  else if (substr(pre,1,2)=="sp") S12 <- matrix.sqrt(var(x))
  Sinv12 <- chol2inv(chol(S12))
    
  if (missing(bgridsize) & binned) bgridsize <- default.gridsize(d)
  if (d > 4) binned <- FALSE
  
  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    ## bin.par <- binning(x.star, bgridsize, sqrt(diag(H.max)))
    ## for large samples, take subset for pilot estimation
    nsub <- n ##min(n, 1e4)
    x.star.sub <- x.star[sample(1:n, size=nsub),]  
    bin.par.sub <- binning(x=x.star.sub, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
  }
  else
    x.star.sub <- x.star
  
  if (d==2)
  {
    if (nstage == 1)
      psi.fun <- psifun1.2d(x.star.sub, pilot=pilot, binned=binned, bin.par=bin.par.sub)
    else if (nstage == 2)
      psi.fun <- psifun2.2d(x.star.sub, pilot=pilot, binned=binned, bin.par=bin.par.sub)
    
    psi40 <- psi.fun[1]
    psi22 <- psi.fun[3]
    psi04 <- psi.fun[5]
    
    ## diagonal bandwidth matrix for 2-dim has exact formula 
    h1 <- (psi04^(3/4)*RK/(psi40^(3/4)*(sqrt(psi40 * psi04)+psi22)*n))^(1/6)
    h2 <- (psi40/psi04)^(1/4) * h1

    return(diag(c(s1^2*h1^2, s2^2*h2^2)))
  }
  else 
  { 
    if (pilot=="amse")
      stop("SAMSE pilot selectors are better for higher dimensions")
 
    ## use normal reference bandwidth as initial condition

    if (missing(Hstart)) 
       Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else    
       Hstart <- Sinv12 %*% Hstart %*% Sinv12

    Hstart <- matrix.sqrt(Hstart)
    
    psi4.mat <- psimat.nd(x.star.sub, d=d, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par.sub)    
  
    ## PI is estimate of AMISE
    pi.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vech(H)) %*% psi4.mat %*% vech(H)
    return(drop(pi.temp)) 
    }
   
    result <- optim(diag(Hstart), pi.temp, method="BFGS")
    H <- diag(result$par) %*% diag(result$par)
  
  ## back-transform
  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
  H <- S12 %*% H %*% S12
  
  return(H)
  }
}



###############################################################################
# Cross-validation bandwidth selectors
###############################################################################

###############################################################################
# Computes the least squares cross validation LSCV function for 2 to 6 dim
# 
# Parameters
# x - data values
# H - bandwidth matrix
#
# Returns
# LSCV(H)
###############################################################################

lscv.1d <- function(x, h, binned=FALSE, bin.par)
{
  n <- length(x)
  lscv1 <- dnorm.sum(x=x, sigma=sqrt(2)*h, inc=1, binned=binned, bin.par=bin.par)
  lscv2 <- dnorm.sum(x=x, sigma=h, inc=0, binned=binned, bin.par=bin.par)
  
  return(lscv1/n^2 - 2/(n*(n-1))*lscv2)     
}

lscv.mat <- function(x, H, binned=FALSE, bin.par)
{
  n <- nrow(x)
  lscv1 <- dmvnorm.sum(x, 2*H, inc=1, binned=binned, bin.par=bin.par)
  lscv2 <- dmvnorm.sum(x, H, inc=0, binned=binned, bin.par=bin.par)
  
  return(lscv1/n^2 - 2/(n*(n-1))*lscv2)     
}

   
###############################################################################
# Finds the bandwidth matrix that minimises LSCV for 2 to 6 dim
# 
# Parameters
# x - data values
# Hstart - initial bandwidth matrix
#
# Returns
# H_LSCV
###############################################################################

hlscv <- function(x, binned=TRUE, bgridsize)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- length(x)
  d <- 1
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))

  if (binned)
  {  
    if (missing(bgridsize)) bgridsize <- 401
    xbin.par <- binning(x, bgridsize=bgridsize, h=hnorm)
    lscv.1d.temp <- function(h)
    {
      return(lscv.1d(x=x, h=h, binned=TRUE, bin.par=xbin.par))
    }
  }
  else
  {
    lscv.1d.temp <- function(h)
    {
      return(lscv.1d(x=x, h=h, binned=FALSE))
    }
  }
  opt <- optimise(f=lscv.1d.temp, interval=c(0.2*hnorm, 5*hnorm, tol=.Machine$double.eps))$minimum
  
  return(opt)
    
}
  
Hlscv <- function(x, Hstart, kfold=1)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")

  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hlscv", kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hlscv", Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  ##RK <- (4*pi)^(-d/2)

  ## use normal reference selector as initial condn
  if (missing(Hstart)) 
    Hstart <- matrix.sqrt((4/ (n*(d + 2)))^(2/(d + 4)) * var(x))
  
  lscv.mat.temp <- function(vechH)
  {
    ##  ensures that H is positive definite
    H <- invvech(vechH) %*% invvech(vechH)
    return(lscv.mat(x=x, H=H, binned=FALSE))
  }
  result <- optim(vech(Hstart), lscv.mat.temp, method="Nelder-Mead")
                                        #control=list(abstol=n^(-10*d)))    
  
  return(invvech(result$par) %*% invvech(result$par))
}

###############################################################################
# Finds the diagonal bandwidth matrix that minimises LSCV for 2 to 6 dim
# 
# Parameters
# x - data values
# Hstart - initial bandwidth matrix
#
# Returns
# H_LSCV,diag
###############################################################################

Hlscv.diag <- function(x, Hstart, binned=FALSE, bgridsize, kfold=1)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")

  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hlscv.diag", kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hlscv.diag", Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  if (missing(Hstart)) 
    Hstart <- matrix.sqrt((4/ (n*(d + 2)))^(2/(d + 4)) * var(x))

  if (missing(bgridsize) & binned) bgridsize <- default.gridsize(d)
  if (d > 4) binned <- FALSE

  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    ## linear binning
    bin.par <- binning(x=x, bgridsize=bgridsize, H=sqrt(diag(diag(H.max))))
  }
  
  lscv.mat.temp <- function(diagH)
  {
    H <- diag(diagH^2)
    return(lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par))
  }
  result <- optim(diag(Hstart), lscv.mat.temp, method="Nelder-Mead")
                 
  return(diag(result$par^2))
}

###############################################################################
# Computes the biased cross validation BCV function for 2-dim
# 
# Parameters
# x - data values
# H1, H2 - bandwidth matrices
#
# Returns
# BCV(H)
###############################################################################

bcv.mat <- function(x, H1, H2)
{
  n <- nrow(x)
  d <- 2

  psi40 <- dmvnorm.deriv.sum(x, Sigma=H2, deriv.order=c(4,0), inc=0)
  psi31 <- dmvnorm.deriv.sum(x, Sigma=H2, deriv.order=c(3,1), inc=0)
  psi22 <- dmvnorm.deriv.sum(x, Sigma=H2, deriv.order=c(2,2), inc=0)
  psi13 <- dmvnorm.deriv.sum(x, Sigma=H2, deriv.order=c(1,3), inc=0)
  psi04 <- dmvnorm.deriv.sum(x, Sigma=H2, deriv.order=c(0,4), inc=0)
    
  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13,psi04)/(n*(n-1))
  psi4.mat <- matrix(coeff * psi.fun, nc=3, nr=3)
  
  RK <- (4*pi)^(-d/2) 
  bcv <- drop(n^(-1)*det(H1)^(-1/2)*RK + 1/4*t(vech(H1)) %*% psi4.mat %*% vech(H1))
  
  return(list(bcv=bcv, psimat=psi4.mat))
}



###############################################################################
# Find the bandwidth matrix that minimises the BCV for 2-dim
# 
# Parameters
# x - data values
# whichbcv - 1 = BCV1
#          - 2 = BCV2 
# Hstart - initial bandwidth matrix
#
# Returns
# H_BCV
###############################################################################

Hbcv <- function(x, whichbcv=1, Hstart, kfold=1)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hbcv", whichbcv=whichbcv, kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hbcv", whichbcv=whichbcv, Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  RK <- (4*pi)^(-d/2)

  # use normal reference b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- Hmax
  ##lo.bound <- -Hmax
  xdiff <- differences(x, upper=FALSE)
  if (missing(Hstart))
    Hstart <- matrix.sqrt(0.9*Hmax)
  bin.par <- binning(x)
  bcv1.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    # ensures that H is positive definite

    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }

  # derivatives of BCV1 function - see thesis  
  bcv1.mat.deriv <- function(vechH)
  {
    H <-  invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    
    psi4.mat <- bcv.mat(x, H, H)$psimat
    psi22 <- psi4.mat[1,3] 
    psi00 <- dmvnorm.sum(x, Sigma=H, inc=0)/(n*(n-1))
    psi22.deriv.xxt <- dmvnorm.deriv.2d.xxt.sum(x, r=c(2,2), Sigma=H)/(n*(n-1))
    psi22.deriv <- t(D2)%*% vec((Hinv %*% psi22.deriv.xxt %*% Hinv + 2* psi00 *Hinv %*% Hinv - psi22*Hinv)/2) 
    
    const <- matrix(c(0,0,1, 0,4,0, 1,0,0), nc=3, byrow=TRUE)
    psi4.mat.deriv<- const %x% psi22.deriv
    
    deriv1 <- -1/2*n^{-1}*RK*t(D2) %*% vec(chol2inv(chol(H)))
    deriv2 <- 1/2 * psi4.mat %*% vech(H) + 1/4 * (t(psi4.mat.deriv) %*% (vech(H) %x% diag(c(1,1,1)))) %*% vech(H)
    
    return(deriv1 + deriv2)      
  }
  
  # derivatives of BCV2 function - see thesis   
  bcv2.mat.deriv <- function(vechH)
  {
    H <-  invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    
    psi4.mat <- bcv.mat(x, H, 2*H)$psimat
    psi22 <- psi4.mat[1,3] 
    psi00 <- dmvnorm.sum(x, Sigma=2*H, inc=0)/(n*(n-1))
    psi22.deriv.xxt <- dmvnorm.deriv.2d.xxt.sum(x,r=c(2,2),Sigma=2*H)/(n*(n-1))
    psi22.deriv <- t(D2)%*% vec((Hinv %*% psi22.deriv.xxt %*% Hinv +
                                 2* psi00 *Hinv %*% Hinv - psi22*Hinv)/2) 
    
    const <- matrix(c(0,0,1, 0,4,0, 1,0,0), nc=3, byrow=TRUE)
    psi4.mat.deriv<- const %x% psi22.deriv
    
    deriv1 <- -1/2*n^{-1}*RK*t(D2) %*% vec(chol2inv(chol(H)))
    deriv2 <- 1/2 * psi4.mat %*% vech(H) + 1/4 *
      (t(psi4.mat.deriv) %*% (vech(H) %x% diag(c(1,1,1)))) %*% vech(H)
    
    return(deriv1 + 2*deriv2)    
  }

  if (whichbcv==1)
    result <- optim(vech(Hstart), bcv1.mat.temp, gr=bcv1.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)))
  else if (whichbcv==2)
    result <- optim(vech(Hstart), bcv2.mat.temp, gr=bcv2.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)))
  
  return(invvech(result$par) %*% invvech(result$par))
}

###############################################################################
# Find the diagonal bandwidth matrix that minimises the BCV for 2-dim
# 
# Parameters
# x - data values
# whichbcv - 1 = BCV1
#          - 2 = BCV2
# Hstart - initial bandwidth matrix
#
# Returns
# H_BCV, diag
###############################################################################

Hbcv.diag <- function(x, whichbcv=1, Hstart, kfold=1)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hbcv.diag", whichbcv=whichbcv, kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hbcv.diag", whichbcv=whichbcv, Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  ##D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  RK <- (4*pi)^(-d/2)
  
  ## use maximally smoothed b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- diag(Hmax)
  ##lo.bound <- rep(0,d)
  
  if (missing(Hstart))
    Hstart <- 0.9*matrix.sqrt(Hmax)

  bcv1.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    ## ensures that H is positive definite

    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }
  
  if (whichbcv == 1)
    result <- optim(diag(Hstart), bcv1.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound))
  else if (whichbcv == 2)
    result <- optim(diag(Hstart), bcv2.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound))

  return(diag(result$par) %*% diag(result$par))
}




########################################################################
### Identifying elements of Theta_6 matrix
########################################################################

Theta6.elem <- function(d)
{
  Theta6.mat <- list()
  Theta6.mat[[d]] <- list()
  for (i in 1:d)
    Theta6.mat[[i]] <- list()
  
  for (i in 1:d)
    for (j in 1:d)
    {  
      temp <- numeric()
      for (k in 1:d)     
        for (ell in 1:d)    
          temp <- rbind(temp, elem(i,d)+2*elem(k,d)+2*elem(ell,d)+elem(j,d))
      
      Theta6.mat[[i]][[j]] <- temp
    }
  
  return(Theta6.mat)
}
  

###############################################################################
# Estimate scalar g_AMSE pilot bandwidth for SCV for 2 to 6 dim
#
# Parameters
# Sigma.star - scaled/ sphered variance matrix
# Hamise - (estimate) of H_AMISE 
# n - sample size
#
# Returns
# g_AMSE pilot bandwidth
###############################################################################


gamse.scv.nd <- function(x.star, d, Sigma.star, Hamise, n, binned=FALSE, bin.par)
{
  psihat6 <- vector()
  g6.star <- gsamse.nd(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(d)

  if (!binned) x.star.diff <- differences(x.star, upper=FALSE)
  ##else bin.par <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(G6.star)))
    
  derivt6 <- deriv.list(d=d, r=6)
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    if (binned)
      psihat6[k] <- kfe(bin.par=bin.par, G=G6.star, r=r, binned=TRUE)
    else 
      psihat6[k] <- kfe.scalar(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }   

  Theta6.mat <- matrix(0, nc=d, nr=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
        temp.sum <- temp.sum + psihat6[which.mat(temp[k,], derivt6)]
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye3 <- diag(d)
  D4 <- dupl(d)$d
  trHamise <- tr(Hamise) ##[1,1] + Hamise[2,2] + Hamise[3,3] 

  ## required constants - see thesis
  Cmu1 <- 1/2*t(D4) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D4)%*% vec(Hamise) + trHamise * t(D4) %*% vec(eye3))

  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) + sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}

###############################################################################
# Estimate unconstrained G_AMSE pilot bandwidth for SCV for 2 to 6 dim
# (J.E. Chacon)
#
# Parameters
# Sigma.star - scaled/ sphered variance matrix
# Hamise - (estimate) of H_AMISE 
# n - sample size
#
# Returns
# G_AMSE pilot bandwidth
###############################################################################

Gamse.scv.nd <- function(x)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  Sd4 <- Sdr(d=d, r=4)  
  Sd6 <- Sdr(d=d, r=6)
  rel.tol<-10^-10  

  ## stage 1 of plug-in
  G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
  ##psi6 <- dmvnorm.deriv.sum(x,Sigma=G6,deriv.order=6,Sdr.mat=Sd6)/n^2
  x.diff <- differences(x, upper=FALSE)
  psihat6 <- kfe(x=x.diff, r=6, G=G6, diff=TRUE)
 
  ## constants for normal reference
  D4phi0 <- drop(dmvnorm.deriv(x=rep(0,d), deriv.order=4, Sdr.mat=Sd4))
  Id1 <- diag(d)
  vId <- vec(Id1)
  Id4 <- diag(d^4)

  ## asymptotic squared bias for r = 4
  AB2r4<-function(vechG){
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=4)%*%D4phi0)*2^(-(d+4)/2)+
      (t(vec(G))%x%Id4)%*%psihat6

    ##return (t(AB) %*% (vec(Hamise) %*% t(vec(Hamise)) %x% diag(d^2)) %*% AB) 
    return (sum(AB^2))
  }

  Hstart <- (4/(d+2))^(2/(d+4))*n^(-2/(d+4))*S
  Hstart <- matrix.sqrt(Hstart)

  res <- optim(vech(Hstart),AB2r4, control=list(reltol=rel.tol),method="BFGS")
  V4 <- res$value
  G4 <- res$par
  G4 <- invvech(G4)%*%invvech(G4)

  return(G4) 
}

###############################################################################
# Computes the smoothed cross validation function for 2 to 6 dim
# 
# Parameters
# x - data values
# H - bandwidth matrix
# G - pilot bandwidth matrix
#
# Returns
# SCV(H)
###############################################################################


scv.1d <- function(x, h, g, binned=TRUE, bin.par, inc=1)
{
  if (!missing(x)) n <- length(x)
  if (!missing(bin.par)) n <- sum(bin.par$counts)
  scv1 <- dnorm.sum(x=x, bin.par=bin.par, sigma=sqrt(2*h^2+2*g^2), binned=binned, inc=inc)
  scv2 <- dnorm.sum(x=x, bin.par=bin.par, sigma=sqrt(h^2+2*g^2), binned=binned, inc=inc)
  scv3 <- dnorm.sum(x=x, bin.par=bin.par, sigma=sqrt(2*g^2), binned=binned, inc=inc)

  bias2 <-  (scv1 - 2*scv2 + scv3)
  if (bias2 < 0) bias2 <- 0
  scv <- (n*h)^(-1)*(4*pi)^(-1/2) + n^(-2)*bias2

  return(scv)
}

scv.mat <- function(x, H, G, binned=FALSE, bin.par, diff=FALSE, n)
{
  if (!diff) n <- nrow(x)  
  d <- ncol(x)

  scv1 <- dmvnorm.sum(x=x, Sigma=2*H + 2*G, inc=1, bin.par=bin.par, binned=binned, diff=diff)
  scv2 <- dmvnorm.sum(x=x, Sigma=H + 2*G, inc=1, bin.par=bin.par, binned=binned, diff=diff)
  scv3 <- dmvnorm.sum(x=x, Sigma=2*G, inc=1, bin.par=bin.par, binned=binned, diff=diff)
  ### need to infer no. x from difference????
  scvmat <- n^(-1)*det(H)^(-1/2)*(4*pi)^(-d/2) + n^(-2)*(scv1 - 2*scv2 + scv3)
    
  return (scvmat)
}



###############################################################################
# Find the bandwidth that minimises the SCV for 1 to 6 dim
# 
# Parameters
# x - data values
# pre - "scale" - pre-scaled data
#     - "sphere"- pre-sphered data
# Hstart - initial bandwidth matrix
#
# Returns
# H_SCV
###############################################################################

hscv <- function(x, nstage=2, binned=TRUE, bgridsize, plot=FALSE)
{
  sigma <- sd(x)
  n <- length(x)
  d <- 1
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))
  if (missing(bgridsize)) bgridsize <- 401
  ##if (missing(hmin))
  hmin <- 0.1*hnorm
  ##if (missing(hmax))
  hmax <- 2*hnorm

  ##bin.par.sub <- binning(x=x[1:min(n, 1e4)], bgridsize=bgridsize, h=hnorm)
  bin.par <- binning(x=x, bgridsize=bgridsize, h=hnorm)
  if (nstage==1)
  {
    psihat6 <- psins.1d(r=6, sigma=sigma)
    psihat10 <- psins.1d(r=10, sigma=sigma)
  }
  else if (nstage==2)
  {
    ##psihat8 <- psins.1d(r=8, sigma=sigma)
    ##psihat12 <- psins.1d(r=12, sigma=sigma)
    g1 <- (2/(7*n))^(1/9)*2^(1/2)*sigma
    g2 <- (2/(11*n))^(1/13)*2^(1/2)*sigma

    psihat6 <- kfe.1d(bin.par=bin.par, binned=TRUE, r=6, g=g1, inc=1)
    psihat10 <- kfe.1d(bin.par=bin.par, binned=TRUE, r=10, g=g2, inc=1)
  }

  g3 <- (-6/((2*pi)^(1/2)*psihat6*n))^(1/7) 
  g4 <- (-210/((2*pi)^(1/2)*psihat10*n))^(1/11)
  psihat4 <- kfe.1d(bin.par=bin.par, binned=TRUE, r=4, g=g3, inc=1)
  psihat8 <- kfe.1d(bin.par=bin.par, binned=TRUE, r=8, g=g4, inc=1)

  C <- (441/(64*pi))^(1/18) * (4*pi)^(-1/5) * psihat4^(-2/5) * psihat8^(-1/9)
  
  scv.1d.temp <- function(h)
  {
    return(scv.1d(x=x, bin.par=bin.par, h=h, g=C*n^(-23/45)*h^(-2), binned=binned, inc=1))
  }

  if (plot)
  {  
    hseq <- seq(hmin,hmax, length=400)
    hscv.seq <- rep(0, length=length(hseq))
    for (i in 1:length(hseq))
      hscv.seq[i] <- scv.1d.temp(hseq[i])
    plot(hseq, hscv.seq, type="l", xlab="h", ylab="SCV(h)")
  }
  
  opt <- optimise(f=scv.1d.temp, interval=c(hmin, hmax))$minimum
  if (n >= 1e5) warning("hscv is not always stable for large samples")
  
  return(opt)
}


Hscv <- function(x, pre="sphere", pilot="samse", Hstart, binned=TRUE, bgridsize, kfold=1)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hscv", pre=pre, pilot=pilot, binned=FALSE, kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hscv", pre=pre, pilot=pilot, binned=FALSE, Hstart=Hstart, kfold=kfold, random=FALSE))
  }

  
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)

  if (substr(pre,1,2)=="sc")
    pre <- "scale"
  else if (substr(pre,1,2)=="sp")
    pre <- "sphere"
 
  if (substr(pilot,1,1)=="a")
    pilot <- "amse"
  else if (substr(pilot,1,1)=="s")
    pilot <- "samse"
  else if (substr(pilot,1,1)=="u")
    pilot <- "unconstr"
  if (pilot=="amse" & d>2)
    stop("SAMSE pilot selectors are better for higher dimensions")
  if (pilot=="unconstr" & d>=6)
    stop("Uconstrained pilots not implemented yet for 6-dim data")

  if(!is.matrix(x)) x <- as.matrix(x)

  ## pre-transform data
  if (pre=="sphere")
    x.star <- pre.sphere(x)
  else if (pre=="scale")
    x.star <- pre.scale(x)
  S.star <- var(x.star)
  n <- nrow(x.star)

  ##if (n > 1000 & !binned)
  ##  warning("Hscv converges slowly for n > 1000 without binned estimation")
  if (d > 4) binned <- FALSE
  if (missing(bgridsize) & binned) bgridsize <- default.gridsize(d)
  
  
  if (pilot=="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x)
    
    Hstart <- matrix.sqrt(Hstart)
    G.amse <- Gamse.scv.nd(x=x)
    x.diff <- differences(x, upper=TRUE)
    
    scv.unconstr.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      scv.temp <- scv.mat(x=x.diff, H=H, G=G.amse, binned=FALSE, diff=TRUE, n=nrow(x))
      return(drop(scv.temp))
    }
    result <- optim(vech(Hstart), scv.unconstr.temp, method="BFGS")
    H <- invvech(result$par) %*% invvech(result$par)
  }
  else if (pilot!="unconstr")
  {
    ##if (binned)
    
    Hamise <- Hpi(x=x, nstage=1, pilot="samse", pre="sphere", binned=binned, bgridsize=bgridsize, kfold=kfold) 
    ##else
    ##  Hamise <- Hpi(x=x, nstage=1, pilot="samse", pre="sphere", binned=FALSE, kfold=kfold)
    
    if (any(is.na(Hamise)))
    {
      warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
      Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    }
    
    if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
    else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
    S12inv <- chol2inv(chol(S12))
    Hamise <- S12inv %*% Hamise%*% S12inv  ## convert to pre-transf data scale
    
    x.star.diff <- differences(x.star, upper=TRUE)
    gamse <- gamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=FALSE)
    G.amse <- gamse^2 * diag(d)

    ## use normal reference bandwidth as initial condition
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else    
      Hstart <- S12inv %*% Hstart %*% S12inv
    Hstart <- matrix.sqrt(Hstart)
    
    scv.mat.temp <- function(vechH)
    {
      H <- invvech(vechH) %*% invvech(vechH)
      return(scv.mat(x.star.diff, H, G.amse, binned=FALSE, diff=TRUE, n=nrow(x)))
    }

    ## back-transform
    result <- optim(vech(Hstart), scv.mat.temp, method= "Nelder-Mead")
    ##control=list(abstol=n^(-10*d)))
    H <- invvech(result$par) %*% invvech(result$par)
    H <- S12 %*% H %*% S12
  }
 
  return(H)
}


Hscv.diag <- function(x, pre="scale", Hstart, binned=FALSE, bgridsize, kfold=1)
{
  if(!is.matrix(x)) x <- as.matrix(x)

  ## k-fold b/w approx
  if (kfold > 1)
  {
    if (missing(Hstart))
      return(Hkfold(x=x, selector="Hscv.diag", pre=pre, binned=FALSE, kfold=kfold, random=FALSE))
    else
      return(Hkfold(x=x, selector="Hscv.diag", pre=pre, binned=FALSE, Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  ## pre-transform data
  
  if (substr(pre,1,2)=="sc")
    pre <- "scale"
  else if (substr(pre,1,2)=="sp")
    pre <- "sphere"

  if (pre=="sphere")
     stop("Using pre-sphering doesn't give a diagonal bandwidth matrix\n")
  
  if (pre=="sphere")
    x.star <- pre.sphere(x)
  else if (pre=="scale")
    x.star <- pre.scale(x)
  
  S.star <- var(x.star)
  n <- nrow(x.star)

  if (missing(bgridsize) & binned) bgridsize <- default.gridsize(d)
  if (d > 4) binned <- FALSE
   
  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))

  S12inv <- chol2inv(chol(S12))
  ##Hamise <- S12inv %*% Hpi.diag(x=x.star, nstage=1, pilot="samse", pre="sphere", binned=binned, bgridsize=bgridsize, kfold=kfold) %*% S12inv
  Hamise <- Hpi.diag(x=x.star, nstage=1, pilot="samse", pre="scale", binned=binned, bgridsize=bgridsize, kfold=kfold)

  if (any(is.na(Hamise)))
  {
    warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
    Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
  }

 
  if (binned)
  {  
    bin.par <- binning(x=x.star, bgridsize=bgridsize, H=diag(diag(Hamise))) 
    gamse <- gamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bin.par=bin.par)
  }
  else
    gamse <- gamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=FALSE)
  G.amse <- gamse^2 * diag(d)
  
  ## use normal reference bandwidth as initial condition
  if (missing(Hstart)) 
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
  else    
    Hstart <- S12inv %*% Hstart %*% S12inv

  Hstart <- matrix.sqrt(Hstart)

  scv.mat.temp <- function(diagH)
  {
    ## ensures that H is positive definite
    H <- diag(diagH) %*% diag(diagH)
    return(scv.mat(x.star, H, G.amse, binned=binned, bin.par=bin.par))
  }
  
  ## back-transform
  result <- optim(diag(Hstart), scv.mat.temp, method= "Nelder-Mead")
  H <- diag(result$par) %*% diag(result$par)
  H <- S12 %*% H %*% S12
  
  return(H)
}



#############################################################################################
## k-fold b/w selectors
#############################################################################################

hkfold <- function(x, selector, k=1, random=FALSE, ...)
{
  n <- length(x)
  d <- 1
  if (random) rand.ind <- sample(1:n)
  else rand.ind <- 1:n

  
  m <- round(n/k,0)
  hstar <- 0
  if (k > 1)
  {
    for (i in 1:(k-1))
    {
      xi <- x[rand.ind[((i-1)*m+1):(i*m)]]
      hi <- do.call(selector, args=list(x=xi, ...))
      hstar <- hstar + hi
    }
  }
  xi <- x[rand.ind[((k-1)*m+1):n]]
  hi <- do.call(selector, args=list(x=xi, ...))
  hstar <- hstar + hi

  hstar <- hstar*k^(-(d+6)/(2*d+8))
  return(hstar)
}

Hkfold <- function(x, selector, kfold=1, random=FALSE, ...)
{
  n <- nrow(x)
  d <- ncol(x)
  if (kfold > n) kfold <- n
  if (random) rand.ind <- sample(1:n)
  else rand.ind <- 1:n

  if (n%%kfold < 10 & n%%kfold >0) kfold <- max(1, kfold-1) 
  m <- round(n/kfold,0)
  Hstar <- matrix(0, ncol=d, nrow=d)
  if (kfold > 1)
  {
    for (i in 1:(kfold-1))
    {
      xi <- x[rand.ind[((i-1)*m+1):(i*m)],]
      Hi <- do.call(selector, args=list(x=xi, ...))
      Hstar <- Hstar + Hi
    }
  }
  xi <- x[rand.ind[((kfold-1)*m+1):n],]
  Hi <- do.call(selector, args=list(x=xi, ...))
  Hstar <- Hstar + Hi

  Hstar <- Hstar*kfold^(-(d+6)/(d+4))
  return(Hstar)
}

