#######################################################################
## Optimal G_r^{NR}
#######################################################################

GNR<-function(r,n,Sigma){
    d<-ncol(Sigma)
    G<-(2/((n*(d+r))))^(2/(d+r+2))*2*Sigma
    return(G)
    }


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
  num <- -2 * dmvnorm.deriv(x=c(0,0), deriv.order=r, Sigma=diag(2), deriv.vec=FALSE)
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

gsamse.nd <- function(Sigma.star, n, modr, nstage=1, psihat=NULL, Sdr.mat)
{
  d <- ncol(Sigma.star)
  K <- numeric(); psi <- numeric()

  ## 4th order g_SAMSE

  K <- dmvnorm.deriv(x=rep(0,d), deriv.order=modr, Sigma=diag(d), add.index=TRUE, deriv.vec=FALSE, Sdr.mat=Sdr.mat)
  K <- K$deriv[apply(K$deriv.ind, 1, is.even)]

  derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE) 
  derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    
  if (modr==4)
  {
    for (i in 1:nrow(derivt4))
    {
      r <- derivt4[i,]
      if (is.even(r))
      {
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
  gamma1 <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  g.samse <- (gamma1 * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}

##############################################################################
## Estimate psi functionals for bivariate data using 1-stage plug-in 
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
##
## Returns
## estimated psi functionals
###############################################################################

psifun1 <- function(x.star, Sd2r4, pilot="samse", binned, bin.par, deriv.order=0, verbose=FALSE, Sdr.mat)
{
  d <- ncol(x.star)
  r <- deriv.order
  S.star <- var(x.star)
  n <- nrow(x.star)
  if (deriv.order>0) pilot <- "vamse"
  
  ## pilots are based on (2r+4)-th order derivatives
  ## compute 1 pilot for VAMSE
  if (pilot=="vamse")
  {
    Sd2r6 <- Sdr(d=d, r=2*r+6)
    D2r4L0 <- DrL0(d=d, r=2*r+4, Sdr=Sd2r4)
    psi2r6.ns <- psins(r=2*r+6, Sigma=S.star, deriv.vec=TRUE, Sdr=Sd2r6) 
    A1 <- sum(D2r4L0^2)
    A2 <- sum(t(D2r4L0) * t(psi2r6.ns) %*% (vec(diag(d)) %x% diag(d^(2*r+4))))
    A3 <- sum((t(psi2r6.ns) %*% (vec(diag(d)) %x% diag(d^(2*r+4))))^2)
    g2r4 <- (((8*r+4*d+16)*A1)/(((-d-2*r-2)*A2 + sqrt((d+2*r+2)^2*A2^2 + (16*r+8*d+32)*A1*A3))*n))^(1/(d+2*r+6))
    psihat.star <- kfe(x=x.star, G=g2r4^2*diag(d), deriv.order=2*r+4, deriv.vec=TRUE, binned=binned, Sdr.mat=Sd2r4, add.index=TRUE, verbose=verbose)
  }  
  ## pilots are based on 4-th order derivatives
  ## compute 1 pilot for SAMSE
  else if (pilot=="samse")
  {
    g.star <- gsamse.nd(S.star, n, 4, Sdr.mat=Sd2r4)
    psihat.star <- kfe(x=x.star, G=g.star^2*diag(d), deriv.order=4, deriv.vec=TRUE, binned=binned, Sdr.mat=Sd2r4,  add.index=TRUE, verbose=verbose)
  }
  ## compute 5 different pilots for AMSE
  else if ((pilot=="amse") & (d==2))
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt4.vec <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=TRUE, only.index=TRUE)

    RK31 <- 15/(64*pi)
    psi00 <- psins(r=c(0,0), Sigma=S.star) 
    psihat.star <- vector()
    g.star <- vector()
    
    for (k in 1:nrow(derivt4))
    { 
      r <- derivt4[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)

      ## odd order
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r=4, n, psi1, psi2, psi00, RK31)
      ## even order
      else
        g.star[k] <- gamse.even.2d(r=4, n, psi1, psi2)[k]
      ##G.star <- g.star[k]^2 * diag(2)

      psihat.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g.star[k], binned=binned, bin.par=bin.par)
    }

    ## create replicated form of psihat
    psihat.star.vec <- rep(0, nrow(derivt4.vec))
    for (k in 1:nrow(derivt4.vec))
      psihat.star.vec[k] <- psihat.star[which.mat(r=derivt4.vec[k,], mat=derivt4)]

    psihat.star <- list(psir=psihat.star.vec, deriv.ind=derivt4.vec)
  }
  
  return(psihat.star)
}


###############################################################################
# Estimate psi functionals for bivariate data using 2-stage plug-in 
#
# Parameters
# x - pre-transformed data points
# pilot - "amse" - different AMSE pilot
#       - "samse" - SAMSE pilot
# Returns
# estimated psi functionals
###############################################################################

psifun2 <- function(x.star, Sd2r4, Sd2r6, pilot="samse", binned, bin.par, deriv.order=0, verbose=FALSE)
{ 
  d <- ncol(x.star)
  r <- deriv.order
  S.star <- var(x.star)
  n <- nrow(x.star)

  ## pilots are based on (2r+4)-th order derivatives
  ## compute 1 pilot for VAMSE
  if (pilot=="vamse")
  {
    Sd2r8 <- Sdr(d=d, r=2*r+8)
    D2r4L0 <- DrL0(d=d, r=2*r+4, Sdr=Sd2r4)
    D2r6L0 <- DrL0(d=d, r=2*r+6, Sdr=Sd2r6)
    psi2r8.ns <- psins(r=2*r+8, Sigma=S.star, deriv.vec=TRUE, Sdr=Sd2r8) 
    A1 <- sum(D2r6L0^2)
    A2 <- sum(t(D2r6L0) * t(psi2r8.ns) %*% (vec(diag(d)) %x% diag(d^(2*r+6))))
    A3 <- sum((t(psi2r8.ns) %*% (vec(diag(d)) %x% diag(d^(2*r+6))))^2)
    g2r6 <- (((8*r+4*d+24)*A1)/(((-d-2*r-4)*A2 + sqrt((d+2*r+4)^2*A2^2 + (16*r+8*d+48)*A1*A3))*n))^(1/(d+2*r+8))
    psi2r6 <- kfe(x=x.star, G=g2r6^2*diag(d), deriv.order=2*r+6, deriv.vec=TRUE, binned=binned, bin.par=bin.par, Sdr.mat=Sd2r6, add.index=FALSE, verbose=verbose)

    A1 <- sum(D2r4L0^2)
    A2 <- sum(t(D2r4L0) * t(psi2r6) %*% (vec(diag(d)) %x% diag(d^(2*r+4))))
    A3 <- sum((t(psi2r6) %*% (vec(diag(d)) %x% diag(d^(2*r+4))))^2)
    g2r4 <- (((8*r+4*d+16)*A1)/(((-d-2*r-2)*A2 + sqrt((d+2*r+2)^2*A2^2 + (16*r+8*d+32)*A1*A3))*n))^(1/(d+2*r+6))
    psihat.star <- kfe(x=x.star, G=g2r4^2*diag(d), deriv.order=2*r+4, deriv.vec=TRUE, binned=binned, bin.par=bin.par, Sdr.mat=Sd2r4, add.index=TRUE, verbose=verbose)
  }
  ## compute 1 pilot for SAMSE    
  else if (pilot=="samse")
  {
    g6.star <- gsamse.nd(S.star, n, 6, Sdr=Sd2r6)
    psihat6.star <- kfe(x=x.star, G=g6.star^2*diag(d), deriv.order=6, deriv.vec=FALSE, binned=binned, bin.par=bin.par, Sdr.mat=Sd2r6, add.index=FALSE, verbose=verbose)
    g.star <- gsamse.nd(S.star, n, 4, nstage=2, psihat=psihat6.star, Sdr=Sd2r4)
    psihat.star <- kfe(x=x.star, G=g.star^2*diag(d), deriv.order=4, deriv.vec=TRUE, binned=binned, bin.par=bin.par, Sdr.mat=Sd2r4, add.index=TRUE, verbose=verbose)
  }
  ## compute different pilots for AMSE
  else if ((pilot=="amse") & (d==2))
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt4.vec <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=TRUE, only.index=TRUE)
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    
    RK31 <- 15/(64*pi)
    RK51 <- 945/(256*pi)
    RK33 <- 225/(256*pi)
    psi00 <- psins(r=rep(0,d), Sigma=S.star) 

    psihat6.star <- vector()
    g6.star <- vector()
    psihat.star <- vector()
    g.star <- vector()
    
    for (k in 1:nrow(derivt6))
    {
      r <- derivt6[k,]
      psi1 <- psins(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins(r=r + 2*elem(2, 2), Sigma=S.star)
      if (prod(r) == 5)
        g6.star[k] <- gamse.odd.2d(r=6, n, psi1, psi2, psi00, RK51)
      else if (prod(r) == 9)
        g6.star[k] <- gamse.odd.2d(r=6, n, psi1, psi2, psi00, RK33) 
      else  
        g6.star[k] <- gamse.even.2d(r=6, n, psi1, psi2)[k]
      
      psihat6.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g6.star[k], binned=binned, bin.par=bin.par) 
    }
    
    ## pilots are based on 4th order derivatives using 6th order psi functionals
    ## computed above 'psihat6.star'
    
    for (k in 1:nrow(derivt4))
    {
      r <- derivt4[k,]
      psi1 <- psihat6.star[7 - (r + 2*elem(1,2))[1]]
      psi2 <- psihat6.star[7 - (r + 2*elem(2,2))[1]]
      
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r=4, n, psi1, psi2, psi00, RK31)
      else
        g.star[k] <- gamse.even.2d(r=4, n, psi1, psi2)[k]

      psihat.star[k] <- kfe.scalar(x=x.star, deriv.order=r, g=g.star[k],  binned=binned, bin.par=bin.par) 
    }

    ## create replicated form of psihat
    psihat.star.vec <- rep(0, nrow(derivt4.vec))
    for (k in 1:nrow(derivt4.vec))
      psihat.star.vec[k] <- psihat.star[which.mat(r=derivt4.vec[k,], mat=derivt4)]

    psihat.star <- list(psir=psihat.star.vec, deriv.ind=derivt4.vec)
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

psifun1.unconstr <- function(x, Sd2r4, binned, bgridsize, deriv.order=0, verbose=FALSE)
{
  n <- nrow(x)
  r <- deriv.order
  S <- var(x)
 
  ## stage 1 of plug-in
  G2r4 <- GNR(r=2*r+4,n=n,Sigma=S) 
  vecPsi2r4 <- kfe(x=x, G=G2r4, deriv.order=2*r+4, binned=binned, bgridsize=bgridsize, double.loop=TRUE, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r4, verbose=verbose) 
  return (vecPsi2r4)
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

psifun2.unconstr <- function(x, Sd2r4, Sd2r6, rel.tol=10^-10, binned, bgridsize, deriv.order=0, verbose=FALSE)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  r <- deriv.order

  ## stage 1 of plug-in
  G2r6 <- GNR(r=2*r+6,n=n,Sigma=S) 
  vecPsi2r6 <- kfe(x=x, G=G2r6, binned=binned, bgridsize=bgridsize, deriv.order=2*r+6, double.loop=TRUE, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r6, verbose=verbose)

  ## asymptotic squared bias for r = 4 for MSE-optimal G
  D2r4phi0 <- DrL0(d=d, r=2*r+4, Sdr=Sd2r4)
  Id2r4 <- diag(d^(2*r+4))
  
  AB2<-function(vechG){
    rr <- 2*r+4
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=rr)%*%D2r4phi0)+(1/2)*(t(vec(G))%x%Id2r4) %*% vecPsi2r6
    return (sum(AB^2))
  }

  Hstart <- GNR(r=2*r+4,n=n,Sigma=S) 
  Hstart <- matrix.sqrt(Hstart)
  res <- optim(vech(Hstart), AB2, control=list(reltol=rel.tol))
  ##V2r4 <- res$value
  G2r4 <- res$par
  G2r4 <- invvech(G2r4)%*%invvech(G2r4)

  ## stage 2 of plug-in
  vecPsi2r4 <- kfe(x=x, G=G2r4, binned=binned, bgridsize=bgridsize, deriv.order=2*r+4, double.loop=TRUE, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r4, verbose=verbose)
  
  return (vecPsi2r4)
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
  
  if (missing(bgridsize)) bgridsize <- default.bgridsize(1)
  return(dpik(x=x, level=nstage, gridsize=bgridsize))
}

Hpi <- function(x, nstage=2, pilot="samse", pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  if (r >0) stop("Currently only deriv.order=0 is implemented")

  Sd2r <- Sdr(d=d,r=2*r)
  vId <- vec(diag(d))
  Idr <- diag(d^r)
  RK <- drop((4*pi)^(-d/2)*2^(-r)*OF(2*r)*Sd2r%*%Kpow(vId,r))

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
    stop("Unconstrained pilots not implemented for >6-dim data")
  
  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  
  Sd2r4 <- Sdr(d=d, r=2*r+4)
  if (nstage==2) Sd2r6 <- Sdr(d=d, r=2*r+6)

  if (pilot=="unconstr")
  {
    ## psi4.mat is on data scale
    ## symmetriser matrices for unconstrained pilot selectors
    if (nstage==1)
      psi.fun <- (-1)^r*psifun1.unconstr(x=x, Sd2r4=Sd2r4, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
    else if (nstage==2)
      psi.fun <- psifun2.unconstr(x=x, Sd2r4=Sd2r4, Sd2r6=Sd2r6, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
    psi2r4.mat <- (-1)^r*invvec(psi.fun)
    
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  else if (pilot!="unconstr")
  {
    if (binned)
    {
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      bin.par.star <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
    }

    ## psi4.mat is on pre-transformed data scale
    if (nstage==1)
      psi.fun <- psifun1(x=x.star, Sd2r4=Sd2r4, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    else if (nstage==2)
      psi.fun <- psifun2(x=x.star, Sd2r4=Sd2r4, Sd2r6=Sd2r6, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    psi2r4.mat <- invvec(psi.fun)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
    else
      Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  ## PI is estimate of AMISE
  pi.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    IdrvH <- Idr%x%vec(H)
    pi.temp <- 1/(det(H)^(1/2)*n)*Kpow(t(vec(Hinv)),r)%*%RK + 1/4* sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))
    return(drop(pi.temp))
  }

  Hstart <- matrix.sqrt(Hstart)
  result <- optim(vech(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
  H <- invvech(result$par) %*% invvech(result$par)
  if (pilot!="unconstr")  H <- S12 %*% H %*% S12     ## back-transform

  if (!amise)
    return(H)
  else
    return(list(H = H, PI.star=result$value))
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


Hpi.diag <- function(x, nstage=2, pilot="samse", pre="scale", Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
{
  if(!is.matrix(x)) x <- as.matrix(x)

  ## k-fold b/w approx
  if (kfold > 1)
  {
    stop("Option kfold > 1 currently disabled in this version")
    ##if (missing(Hstart))
    ##  return(Hkfold(x=x, selector="Hpi.diag", kfold=kfold, random=FALSE, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE, amise=amise))
    ##else
    ##  return(Hkfold(x=x, selector="Hpi.diag", kfold=kfold, random=FALSE, Hstart=Hstart, nstage=nstage, pilot=pilot, pre=pre, binned=FALSE, amise=amise))
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
  r <- deriv.order
  if (r >0) stop("Currently only deriv.order=0 is implemented")
  
  if (substr(pre,1,2)=="sc") S12 <- diag(sqrt(diag(var(x))))
  else if (substr(pre,1,2)=="sp") S12 <- matrix.sqrt(var(x))
  Sinv12 <- chol2inv(chol(S12))

  if (d > 4) binned <- FALSE
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par <- binning(x=x.star, bgridsize=bgridsize, H=diag(H.max)) 
  }
 
  Sd2r4 <- Sdr(d=d, r=2*r+4)
  if (nstage==2) Sd2r6 <- Sdr(d=d, r=2*r+6)

  if (d==2)
  {
    if (nstage == 1)
      psi.fun <- psifun1(x.star, pilot=pilot, binned=binned, bin.par=bin.par, verbose=verbose, deriv.order=r, Sd2r4=Sd2r4)$psir
    else if (nstage == 2)
      psi.fun <- psifun2(x.star, pilot=pilot, binned=binned, bin.par=bin.par, verbose=verbose, deriv.order=r, Sd2r4=Sd2r4, Sd2r6=Sd2r6)$psir

    psi40 <- psi.fun[1]
    psi22 <- psi.fun[6]
    psi04 <- psi.fun[16]
    
    ## diagonal bandwidth matrix for 2-dim has exact formula 
    h1 <- (psi04^(3/4)*RK/(psi40^(3/4)*(sqrt(psi40*psi04)+psi22)*n))^(1/6)
    h2 <- (psi40/psi04)^(1/4) * h1

    H <- diag(c(s1^2*h1^2, s2^2*h2^2))

    psimat4.D <- invvech(c(psi40, psi22, psi04))
    amise.star <- drop(n^(-1)*RK*(h1*h2)^(-1) + 1/4*c(h1,h2)^2 %*% psimat4.D %*% c(h1,h2)^2)
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

    if (nstage == 1)
      psi.fun <- psifun1(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose, Sd2r4=Sd2r4)$psir
    else if (nstage == 2)
      psi.fun <- psifun2(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose, Sd2r4=Sd2r4, Sd2r6=Sd2r6)$psir
    psi4.mat <- invvec(psi.fun)
  
    ## PI is estimate of AMISE
    pi.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vec(H)) %*% psi4.mat %*% vec(H)
    return(drop(pi.temp)) 
    }

    result <- optim(diag(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
  
    ## back-transform
    if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
    else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
    amise.star <- result$value
    H <- S12 %*% H %*% S12
  }

  if (!amise)
    return(H)
  else
    return(list(H = H, PI.star=amise.star))
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

lscv.1d <- function(x, h, binned, bin.par, deriv.order=0)
{
  r <- deriv.order
 
  lscv1 <- kfe.1d(x=x, g=sqrt(2)*h, inc=1, binned=binned, bin.par=bin.par, deriv.order=2*r) 
  lscv2 <- kfe.1d(x=x, g=h, inc=0, binned=binned, bin.par=bin.par, deriv.order=2*r) 
  return((-1)^r*(lscv1 - 2*lscv2))     
}

lscv.mat <- function(x, H, binned=FALSE, bin.par, bgridsize, deriv.order=0, verbose=FALSE, Sd2r, kfold)
{
  r <- deriv.order
  d <- ncol(x)
  
  lscv1 <- kfe(x=x, G=2*H, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, Sdr.mat=Sd2r, verbose=verbose, kfold=kfold)$psir
  lscv2 <- kfe(x=x, G=H, inc=0, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, Sdr.mat=Sd2r, verbose=verbose, kfold=kfold)$psir
  lscv <- drop(lscv1 - 2*lscv2)
  lscv <- (-1)^r*sum(vec(diag(d^r))*lscv)  
  return(lscv)  
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

hlscv <- function(x, binned=TRUE, bgridsize, deriv.order=0)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- length(x)
  d <- 1
  r <- deriv.order
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))

  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (binned) bin.par <- binning(x, bgridsize=bgridsize, h=hnorm)
  lscv.1d.temp <- function(h)
  {
    return(lscv.1d(x=x, h=h, binned=binned, bin.par=bin.par, deriv.order=r))
  }
  opt <- optimise(f=lscv.1d.temp, interval=c(0.2*hnorm, 5*hnorm, tol=.Machine$double.eps))$minimum
  
  return(opt)
    
}
  
Hlscv <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")

  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order 

  ## use normal reference selector as initial condn
  if (missing(Hstart)) Hstart <- matrix.sqrt((4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  
  if (d > 4) binned <- FALSE
  Sd2r <- Sdr(d=d, r=2*r)
  
  lscv.mat.temp <- function(vechH)
  {
    ##  ensures that H is positive definite
    H <- invvech(vechH) %*% invvech(vechH)
    val <- lscv.mat(x=x, H=H, binned=binned, bgridsize=bgridsize, deriv.order=r, Sd2r=Sd2r, verbose=FALSE, kfold=kfold)
    return(val)
  }
  result <- optim(vech(Hstart), lscv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))     

  H <- invvech(result$par) %*% invvech(result$par)
  amise.opt <- result$value

  if (!amise)
    return(H)
  else
    return(list(H=H, LSCV=amise.opt))
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

Hlscv.diag <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")

  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  
  if (missing(Hstart)) 
    Hstart <- matrix.sqrt((4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x))
  
  if (d > 4) binned <- FALSE

  ## linear binning
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    RK <- (4*pi)^(-d/2)
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    bin.par <- binning(x=x, bgridsize=bgridsize, H=sqrt(diag(diag(H.max))))
  }
  
  lscv.mat.temp <- function(diagH)
  {
    H <- diag(diagH^2)
    return(lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par, deriv.order=r, verbose=FALSE, kfold=kfold))
  }
  result <- optim(diag(Hstart), lscv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
  H <- diag(result$par^2)
  amise.opt <- result$value

  if (!amise)
    return(H)
  else
    return(list(H=H, LSCV=amise.opt))
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

  psi <- kfe(x, G=H2, deriv.order=4, add.index=TRUE, deriv.vec=TRUE, inc=0)
  psi40 <- psi$psir[1]
  psi31 <- psi$psir[2]
  psi22 <- psi$psir[4]
  psi13 <- psi$psir[8]
  psi04 <- psi$psir[16]
  
  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13,psi04)/(n*(n-1))
  psi4.mat <- matrix(coeff * psi.fun, ncol=3, nrow=3)
  
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

Hbcv <- function(x, whichbcv=1, Hstart, amise=FALSE, kfold=1, verbose=FALSE)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    stop("Option kfold > 1 currently disabled in this version")
    ##if (missing(Hstart))
    ##  return(Hkfold(x=x, selector="Hbcv", whichbcv=whichbcv, kfold=kfold, random=FALSE))
    ##else
    ##  return(Hkfold(x=x, selector="Hbcv", whichbcv=whichbcv, Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  ##D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  RK <- (4*pi)^(-d/2)

  ## use normal reference b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- Hmax
  if (missing(Hstart)) Hstart <- matrix.sqrt(0.9*Hmax)

  bcv1.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    ## ensures that H is positive definite

    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }

  if (whichbcv==1)
    result <- optim(vech(Hstart), bcv1.mat.temp, #gr=bcv1.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)), control=list(trace=as.numeric(verbose)))
  else if (whichbcv==2)
    result <- optim(vech(Hstart), bcv2.mat.temp, #gr=bcv2.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)), control=list(trace=as.numeric(verbose)))
  H <- invvech(result$par) %*% invvech(result$par)
  amise.opt <- result$value

  if (!amise)
    return(H)
  else
    return(list(H = H, BCV=amise.opt))
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

Hbcv.diag <- function(x, whichbcv=1, Hstart, amise=FALSE, kfold=1, verbose=FALSE)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    stop("Option kfold > 1 currently disabled in this version")
    ##if (missing(Hstart))
    ##  return(Hkfold(x=x, selector="Hbcv.diag", whichbcv=whichbcv, kfold=kfold, random=FALSE))
    ##else
    ##  return(Hkfold(x=x, selector="Hbcv.diag", whichbcv=whichbcv, Hstart=Hstart, kfold=kfold, random=FALSE))
  }
  
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  ## use maximally smoothed b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- diag(Hmax)
  
  if (missing(Hstart))
    Hstart <- 0.9*matrix.sqrt(Hmax)

  bcv1.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }
  
  if (whichbcv == 1)
    result <- optim(diag(Hstart), bcv1.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound), control=list(trace=as.numeric(verbose)))
  else if (whichbcv == 2)
    result <- optim(diag(Hstart), bcv2.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound), control=list(trace=as.numeric(verbose)))

  H <- diag(result$par) %*% diag(result$par)
  amise.opt <- result$value
  if (!amise)
    return(H)
  else
    return(list(H = H, BCV=amise.opt))
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

gamse.scv.nd <- function(x.star, Sd6, d, Sigma.star, Hamise, n, binned=FALSE, bin.par, bgridsize, verbose=FALSE, nstage=1)
{
  if (nstage==0)
  {
    psihat6.star <- psins(r=6, Sigma=Sigma.star, deriv.vec=TRUE, Sdr=Sd6) 
  }
  else if (nstage==1)
  {  
    g6.star <- gsamse.nd(Sigma.star, n, 6) 
    G6.star <- g6.star^2 * diag(d)
    psihat6.star <- kfe(x=x.star, bin.par=bin.par, Sdr.mat=Sd6, deriv.order=6, G=G6.star, deriv.vec=FALSE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }

  derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
  Theta6.mat <- matrix(0, ncol=d, nrow=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
        temp.sum <- temp.sum + psihat6.star[which.mat(temp[k,], derivt6)]
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye3 <- diag(d)
  D4 <- dupl(d)$d
  trHamise <- tr(Hamise) 

  ## required constants - see thesis
  Cmu1 <- 1/2*t(D4) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D4)%*% vec(Hamise) + trHamise * t(D4) %*% vec(eye3))

  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) + sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}

gvamse.scv.nd <- function(x.star, Sd6, Sd4, d, Sigma.star, Hamise, n, binned=FALSE, bin.par, bgridsize, verbose=FALSE, nstage=1)
{
  if (nstage==0)
  {
    psihat6.star <- psins(r=6, Sigma=Sigma.star, deriv.vec=TRUE, Sdr=Sd6) 
  }
  else if (nstage==1)
  {  
    g6.star <- gsamse.nd(Sigma.star, n, 6) 
    G6.star <- g6.star^2 * diag(d)
    psihat6.star <- kfe(x=x.star, bin.par=bin.par, Sdr.mat=Sd6, deriv.order=6, G=G6.star, deriv.vec=TRUE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }
  
  ## constants for normal reference
  D4phibar0 <- drop(dmvnorm.deriv(x=rep(0,d), Sigma=2*diag(d), deriv.order=4, Sdr.mat=Sd4))
  Id1 <- diag(d)
  vId <- vec(Id1)
  Id2 <- diag(d^2)
  Id4 <- diag(d^4)

  A1 <- t(D4phibar0) %*% ((vec(Hamise) %*% t(vec(Hamise))) %x% Id2) %*% D4phibar0
  A2 <- t(D4phibar0) %*% ((vec(Hamise) %*% t(vec(Hamise))) %x% Id2) %*% (t(vId) %x% Id4) %*% psihat6.star
  A3 <- t(psihat6.star) %*% (vId %x% Id4) %*% ((vec(Hamise) %*% t(vec(Hamise))) %x% Id2) %*% (t(vId) %x% Id4) %*% psihat6.star

  gvamse <- (2*(d+4)*A1/((-(d+2)*A2 + sqrt((d+2)^2*A2^2 + 8*(d+4)*A1*A3))*n))^(1/(d+6))
  return(drop(gvamse))
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

Gamse.scv.nd <- function(x, Sd6, Sd4, binned=FALSE, bin.par, bgridsize, rel.tol=10^-10, verbose=FALSE, nstage=1)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  ## stage 1 of plug-in
  if (nstage==1)
  {  
    G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
    psihat6 <- kfe(x=x, deriv.order=6, Sdr.mat=Sd6, G=G6, deriv.vec=TRUE, add.index=FALSE, binned=binned, bin.par=bin.par, bgridsize=bgridsize, verbose=verbose)
  }
  else if (nstage==0)
  {
    psihat6 <- psins(r=6, Sigma=S, deriv.vec=TRUE, Sdr=Sd6) 
  }
  
  ## constants for normal reference
  D4phi0 <- drop(dmvnorm.deriv(x=rep(0,d), deriv.order=4, Sdr.mat=Sd4))
  ##Id1 <- diag(d)
  ##vId <- vec(Id1)
  Id4 <- diag(d^4)

  ## asymptotic squared bias for r = 4
  AB2r4<-function(vechG){
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=4)%*%D4phi0)*2^(-(d+4)/2) + (t(vec(G))%x%Id4)%*%psihat6
    return (sum(AB^2))
  }

  Hstart <- (4/(d+2))^(2/(d+4))*n^(-2/(d+4))*S
  Hstart <- matrix.sqrt(Hstart)

  res <- optim(vech(Hstart), AB2r4, control=list(reltol=rel.tol, trace=as.numeric(verbose)), method="BFGS")
  ##V4 <- res$value
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


scv.1d <- function(x, h, g, binned=TRUE, bin.par, inc=1, deriv.order=0)
{
  r <- deriv.order
  if (!missing(x)) n <- length(x)
  if (!missing(bin.par)) n <- sum(bin.par$counts)
  scv1 <- kfe.1d(x=x, deriv.order=r, bin.par=bin.par, g=sqrt(2*h^2+2*g^2), binned=binned, inc=inc)
  scv2 <- kfe.1d(x=x, deriv.order=r, bin.par=bin.par, g=sqrt(h^2+2*g^2), binned=binned, inc=inc)
  scv3 <- kfe.1d(x=x, deriv.order=r, bin.par=bin.par, g=sqrt(2*g^2), binned=binned, inc=inc)

  bias2 <-  (scv1 - 2*scv2 + scv3)
  if (bias2 < 0) bias2 <- 0
  scv <- (n*h)^(-1)*(4*pi)^(-1/2) + bias2

  return(scv)
}

scv.mat <- function(x, H, G, binned=FALSE, bin.par, bgridsize, verbose=FALSE, deriv.order=0, Sdr.mat)
{
  d <- ncol(x)
  n <- nrow(x)
  r <- deriv.order
  scv1 <- kfe(x=x, G=2*H + 2*G, deriv.order=r, Sdr.mat=Sdr.mat, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose)$psir
  scv2 <- kfe(x=x, G=H + 2*G, deriv.order=r, Sdr.mat=Sdr.mat, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose)$psir
  scv3 <- kfe(x=x, G=2*G, deriv.order=r, Sdr.mat=Sdr.mat, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose)$psir

  bias2 <- scv1 - 2*scv2 + scv3
  if (bias2 < 0) bias2 <- 0
  scvmat <- n^(-1)*det(H)^(-1/2)*(4*pi)^(-d/2) + bias2
     
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
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  hmin <- 0.1*hnorm
  hmax <- 2*hnorm

  bin.par <- binning(x=x, bgridsize=bgridsize, h=hnorm)
  if (nstage==1)
  {
    psihat6 <- psins.1d(r=6, sigma=sigma)
    psihat10 <- psins.1d(r=10, sigma=sigma)
  }
  else if (nstage==2)
  {
    g1 <- (2/(7*n))^(1/9)*2^(1/2)*sigma
    g2 <- (2/(11*n))^(1/13)*2^(1/2)*sigma

    psihat6 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=6, g=g1, inc=1)
    psihat10 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=10, g=g2, inc=1)
  }

  g3 <- (-6/((2*pi)^(1/2)*psihat6*n))^(1/7) 
  g4 <- (-210/((2*pi)^(1/2)*psihat10*n))^(1/11)
  psihat4 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=4, g=g3, inc=1)
  psihat8 <- kfe.1d(x=x, bin.par=bin.par, binned=binned, deriv.order=8, g=g4, inc=1)

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


Hscv <- function(x, nstage=2, pre="sphere", pilot="samse", Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, verbose=FALSE)
{
  ## k-fold b/w approx
  if (kfold > 1)
  {
    stop("Option kfold > 1 currently disabled in this version")
    ##if (missing(Hstart))
    ##  return(Hkfold(x=x, selector="Hscv", pre=pre, pilot=pilot, binned=FALSE, kfold=kfold, random=FALSE))
    ##else
    ##  return(Hkfold(x=x, selector="Hscv", pre=pre, pilot=pilot, binned=FALSE, Hstart=Hstart, kfold=kfold, random=FALSE))
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
  else if (substr(pilot,1,1)=="v")
    pilot <- "vamse"
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

  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)

  Sd0 <- Sdr(d=d, r=0)
  Sd4 <- Sdr(d=d, r=4)
  Sd6 <- Sdr(d=d, r=6)

  if (pilot=="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x)

    Hstart <- matrix.sqrt(Hstart)
    G.amse <- Gamse.scv.nd(x=x, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, Sd4=Sd4, nstage=nstage-1) 

    scv.unconstr.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      scv.temp <- scv.mat(x=x, H=H, G=G.amse, binned=binned, bgridsize=bgridsize, verbose=FALSE, Sdr.mat=Sd0)
      return(drop(scv.temp))
    }
    result <- optim(vech(Hstart), scv.unconstr.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
  }
  else if (pilot!="unconstr")
  {
    Hamise <- Hpi(x=x, nstage=1, pilot=pilot, pre="sphere", binned=binned, bgridsize=bgridsize, verbose=verbose) 
    if (any(is.na(Hamise)))
    {
      warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
      Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    }
  
    if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
    else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
    S12inv <- chol2inv(chol(S12))
    Hamise <- S12inv %*% Hamise%*% S12inv  ## convert to pre-transf data scale

    if (pilot=="amse" | pilot=="samse")
      gamse <- gamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, nstage=nstage-1)
    else if (pilot=="vamse")
      gamse <- gvamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, Sd4=Sd4, nstage=nstage-1)
    G.amse <- gamse^2 * diag(d)

    scv.mat.temp <- function(vechH)
    {
      H <- invvech(vechH) %*% invvech(vechH)
      return(scv.mat(x.star, H, G.amse, binned=binned, bgridsize=bgridsize, verbose=FALSE, Sdr.mat=Sd0))
    }

    ## use normal reference bandwidth as initial condition
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else    
      Hstart <- S12inv %*% Hstart %*% S12inv
    Hstart <- matrix.sqrt(Hstart)

    result <- optim(vech(Hstart), scv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    H <- S12 %*% H %*% S12  ## back-transform
  }

  if (!amise)
    return(H)
  else
    return(list(H = H, SCV.star=result$value))
}


Hscv.diag <- function(x, nstage=2, pre="scale", pilot="samse", Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, verbose=FALSE)
{
  if(!is.matrix(x)) x <- as.matrix(x)

  ## k-fold b/w approx
  if (kfold > 1)
  {
    stop("Option kfold > 1 currently disabled in this version")
    ##if (missing(Hstart))
    ##  return(Hkfold(x=x, selector="Hscv.diag", pre=pre, binned=FALSE, kfold=kfold, random=FALSE))
    ##else
    ##  return(Hkfold(x=x, selector="Hscv.diag", pre=pre, binned=FALSE, Hstart=Hstart, kfold=kfold, random=FALSE))
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

  if (d > 4) binned <- FALSE
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
  }  

  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))

  S12inv <- chol2inv(chol(S12))
  Hamise <- Hpi.diag(x=x.star, nstage=1, pilot="samse", pre="scale", binned=binned, bgridsize=bgridsize, verbose=verbose)
  ##if (verbose) cat("Plugin stage complete.\n")
  
  if (any(is.na(Hamise)))
  {
    warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
    Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
  }

  Sd4 <- Sdr(d=d, r=4)
  Sd6 <- Sdr(d=d, r=6)

  if (pilot=="amse" | pilot=="samse")
    gamse <- gamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, nstage=nstage-1)
  else if (pilot=="vamse")
    gamse <- gvamse.scv.nd(x.star=x.star, d=d, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, Sd4=Sd4, nstage=nstage-1)
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
    return(scv.mat(x.star, H, G.amse, binned=binned, bin.par=bin.par, verbose=FALSE))
  }
  
  ## back-transform
  result <- optim(diag(Hstart), scv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
  H <- diag(result$par) %*% diag(result$par)
  H <- S12 %*% H %*% S12

  if (!amise)
    return(H)
  else
    return(list(H = H, SCV.star=result$value))
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




