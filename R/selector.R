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

gsamse <- function(Sigma.star, n, modr, nstage=1, psihat=NULL, Sdr.mat)
{
  d <- ncol(Sigma.star)
  K <- numeric(); psi <- numeric()

  ## 4th order g_SAMSE

  K <- dmvnorm.deriv(x=rep(0,d), deriv.order=modr, Sigma=diag(d), add.index=TRUE, deriv.vec=FALSE, Sdr.mat=Sdr.mat)
  K <- K$deriv[apply(K$deriv.ind, 1, is.even)]
 
  if (modr==4)
  {
    derivt4 <- dmvnorm.deriv(x=rep(0,d), deriv.order=4, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
  
    for (i in 1:nrow(derivt4))
    {
      r <- derivt4[i,]
      if (is.even(r))
      {
        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
          {
            A3psi <- A3psi + psins(r=r+2*elem(j,d), Sigma=Sigma.star)
            ##A3psi <- A3psi + psins6[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          }
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
    derivt6 <- dmvnorm.deriv(x=rep(0,d), deriv.order=6, add.index=TRUE, deriv.vec=FALSE, only.index=TRUE)
  
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



## Scalar pilot selector for derivatives r>0 from Chacon & Duong (2011)
gdscalar <- function(x, d, r, n, verbose, nstage=1, scv=FALSE)
{
  if (scv) cf <- c(2^(-d), 2^(-d/2+1), 4)
  else cf <- c(1,1,1)
  
  if (nstage==1)
  {
    G2r4 <- GNR(r=2*r+4,n=n,Sigma=var(x))
    g2r4 <- sqrt(G2r4[1,1])
  }
  else if (nstage==2)
  {
    G2r6.NR <- GNR(r=2*r+6,n=n,Sigma=var(x))
    g2r6.nr <- prod(sqrt(diag(G2r6.NR)))^(1/d)##sqrt(G2r6.NR[1,1])
    L0 <- dmvnorm.mixt(x=rep(0,d), mus=rep(0,d), Sigmas=diag(d), props=1)
    eta2r6 <- eta.kfe.y(x=x, deriv.order=2*r+6, G=g2r6.nr^2*diag(d), verbose=verbose, symm=FALSE)
    A1 <- cf[1]*(2*d+4*r+8)*L0^2*OF(2*r+4)*nu(r=r+2, A=diag(d))
    A2 <- cf[2]*(-1)^(r+2)*(d+2*r+2)*L0*OF(2*r+4)*eta2r6
    A3 <- cf[3]*eta2r6^2
    
    g2r4 <- (2*A1/((-A2+ sqrt(A2^2 +4*A1*A3))*n))^(1/(d+2*r+6))
  }
  return(g2r4)
}

##############################################################################
## Scalar pilot selector for derivatives r>0 from Chacon & Duong (2011)
## Generalisation of gsamse for r>0
###############################################################################


## Unconstrained pilot selector for derivatives r>0 from Chacon & Duong (2011)
## Generalisation of Gunconstr for r>0

Gdunconstr <- function(x, d, r, n, nstage=1, verbose, scv=FALSE, Sdr.flag=TRUE, optim.fun="nlm", thin=1)
{
  if (scv) cf <- c(2^(-d/2), 2)
  else cf <- c(1,1)
  ##if (scv) cf <- c(2^(-d), 2^(-d/2+1), 4)
  ##else cf <- c(1,1,1)
  S <- var(x)

  if (nstage==1)
  {
    G2r4 <- GNR(r=2*r+4,n=n,Sigma=S)
  }
  else if (nstage==2)
  {
    G2r4.NR <- GNR(r=2*r+4,n=n,Sigma=S)
    G2r6.NR <- GNR(r=2*r+6,n=n,Sigma=S)
    vecPsi2r6 <- kfe(x=x, G=G2r6.NR, binned=FALSE, deriv.order=2*r+6, deriv.vec=TRUE, add.index=FALSE, verbose=verbose, Sdr.flag=Sdr.flag, thin=thin)   
    D2r4phi0 <- DrL0(d=d, r=2*r+4, Sdr.flag=Sdr.flag, verbose=verbose, thin=thin)
    ##Id2r4 <- diag(d^(2*r+4))
    dls <- (0:(d^2-1))*d^(2*r+4) 

    AB2 <- function(vechG)
    {
      G <- invvech(vechG) %*% invvech(vechG)
      Ginv <- chol2inv(chol(G))
      Ginv12 <- matrix.sqrt(Ginv)
      
      ###eta functionals are more efficient but are currently less accurate than direct computation of psi functionals (kfe)
      ##bias1 <- cf[1]*n^(-2)*det(G)^(-1)*L0^2*OF(2*r+4)*nu(r=r+2, A=Ginv%*%Ginv)
      ##bias2 <- cf[2]*(-1)^(r+2)*n^(-1)*det(G)^(-1/2)*L0*OF(2*r+4)*eta.rs.kfe(x=x, r=2, s=2*r+4, A=G, B=Ginv, C=G2r6.NR)    
      ##bias3 <- cf[3]*1/4*eta.rs.kfe(x=x, r=2, s=2*r+4, A=G, B=diag(d), C=G2r6.NR)^2
      ##bias3 <- cf[3]*1/4*eta2r6*eta.rs.kfe(x=x, r=4, s=2*r+2, A=G, B=diag(d), C=G2r6.NR)
      ##AB2.val <- bias1 + bias2 + bias3
      
      ## direct computation
      v1 <- n^(-1)*det(Ginv12)*drop(Kpow(A=Ginv12,pow=2*r+4)%*%D2r4phi0) 
      ##v2 <- (1/2)*drop((t(vec(G))%x%Id2r4) %*% vecPsi2r6)
      
      ## direct, recursive computation by Jose Chacon 02/2012
      ##v1 <- n^(-1)*dmvnorm(rep(0,d),rep(0,d),sigma=G)*(-1)^(r+2)*OF(2*r+4)*Sdrv(d=d,r=2*r+4,v=Kpow(vec(Ginv),r+2))
      v2 <- numeric(d^(2*r+4))
        for(k in 1:d^(2*r+4)){v2[k]<-(1/2)*sum(vec(G)*vecPsi2r6[dls+k])}
      
      AB <- cf[1]*v1 + cf[2]*v2
      AB2.val <- sum(AB^2)
      return(AB2.val)
    }
    Gstart <- matrix.sqrt(G2r4.NR)
    optim.fun1 <- tolower(substr(optim.fun, 1,1))
    if (optim.fun1=="n")
    {
      result <- nlm(p=vech(Gstart), f=AB2, print.level=2*as.logical(verbose))
      G2r4 <- result$estimate
    } 
    else
    { 
      result <- optim(vech(matrix.sqrt(Gstart)), AB2, method="BFGS", control=list(trace=as.numeric(verbose)))
      G2r4 <- result$par
    }   
    G2r4 <- invvech(G2r4)%*%invvech(G2r4)
  }
  return(G2r4)
}

##############################################################################
## Estimate psi functionals using 1-stage plug-in 
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
##
## Returns
## estimated psi functionals
###############################################################################

psifun1 <- function(x.star, Sd2r4, pilot="samse", binned, bin.par, deriv.order=0, verbose=FALSE)
{
  d <- ncol(x.star)
  r <- deriv.order
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  ## pilots are based on (2r+4)-th order derivatives
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
  {
    g.star <- gsamse(S.star, n, 4, Sdr.mat=Sd2r4)
    psihat.star <- kfe(x=x.star, G=g.star^2*diag(d), deriv.order=4, deriv.vec=TRUE, binned=binned, Sdr.mat=Sd2r4, add.index=TRUE, verbose=verbose)
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
      if (prod(r) == 3)g.star[k] <- gamse.odd.2d(r=4, n, psi1, psi2, psi00, RK31)
      ## even order
      else  g.star[k] <- gamse.even.2d(r=4, n, psi1, psi2)[k]
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
# Estimate psi functionals using 2-stage plug-in 
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
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
  {
    g6.star <- gsamse(S.star, n=n, modr=6, Sdr.mat=Sd2r6)
    psihat6.star <- kfe(x=x.star, G=g6.star^2*diag(d), deriv.order=6, deriv.vec=TRUE, binned=binned, bin.par=bin.par, Sdr.mat=Sd2r6, add.index=FALSE, verbose=verbose)
    g.star <- gsamse(S.star, n=n, modr=4, nstage=2, psihat=psihat6.star)
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

psifun1.unconstr <- function(x, Sd2r4, binned, bgridsize, deriv.order=0, verbose=FALSE, Sdr.flag=TRUE)
{
  n <- nrow(x)
  r <- deriv.order
  S <- var(x)
 
  ## stage 1 of plug-in
  G2r4 <- GNR(r=2*r+4,n=n,Sigma=S) 

  vecPsi2r4 <- kfe(x=x, G=G2r4, deriv.order=2*r+4, binned=binned, bgridsize=bgridsize, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r4, verbose=verbose, Sdr.flag=Sdr.flag) 
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

psifun2.unconstr <- function(x, Sd2r4, Sd2r6, rel.tol=10^-10, binned, bgridsize, deriv.order=0, verbose=FALSE, Sdr.flag=TRUE, optim.fun="nlm")
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  r <- deriv.order

  ## stage 1 of plug-in
  G2r6 <- GNR(r=2*r+6,n=n,Sigma=S) 
  vecPsi2r6 <- kfe(x=x, G=G2r6, binned=binned, bgridsize=bgridsize, deriv.order=2*r+6, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r6, verbose=verbose, Sdr.flag=Sdr.flag)

  ## asymptotic squared bias for r = 4 for MSE-optimal G
  D2r4phi0 <- DrL0(d=d, r=2*r+4, Sdr.flag=Sdr.flag)
  Id2r4 <- diag(d^(2*r+4))
  
  AB2<-function(vechG){
    rr <- 2*r+4
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=rr)%*%D2r4phi0)+(1/2)*(t(vec(G))%x%Id2r4) %*% vecPsi2r6
    return (sum(AB^2))
  }
      
  Gstart <- GNR(r=2*r+4,n=n,Sigma=S) 
  Gstart <- matrix.sqrt(Gstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    res <- nlm(p=vech(Gstart), f=AB2, print.level=2*as.logical(verbose))    
    G2r4 <- res$estimate
  }
  else 
  {
    res <- optim(vech(Gstart), AB2, control=list(reltol=rel.tol, trace=as.numeric(verbose)))
    G2r4 <- res$par
  }
  G2r4 <- invvech(G2r4)%*%invvech(G2r4)

  ## stage 2 of plug-in
  vecPsi2r4 <- kfe(x=x, G=G2r4, binned=binned, bgridsize=bgridsize, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, Sdr.mat=Sd2r4, verbose=verbose, Sdr.flag=Sdr.flag)
  
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

Hpi <- function(x, nstage=2, pilot="samse", pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", Sdr.flag=TRUE)
{
  n <- nrow(x)
  d <- ncol(x)
  ##S <- var(x)
  r <- deriv.order

  if(!is.matrix(x)) x <- as.matrix(x)
  if (substr(pilot,1,1)=="a") pilot <- "amse"
  else if (substr(pilot,1,1)=="s") pilot <- "samse"
  else if (substr(pilot,1,1)=="u") pilot <- "unconstr"
  else if (substr(pilot,1,2)=="du") pilot <- "dunconstr"                     
  else if (substr(pilot,1,2)=="ds") pilot <- "dscalar"
 
  if (pilot=="amse" & (d>2 | r>0)) stop("AMSE pilot selectors not defined for d>2 and/or r>0.")
  if ((pilot=="samse" | pilot=="unconstr") & r>0) stop("DSCALAR or DUNCONSTR pilot selectors are better for derivatives r>0.")
  if (pilot=="unconstr" & d>=6) stop("Unconstrained pilots are not implemented for d>6.")
  
  if (substr(pre,1,2)=="sc") pre <- "scale"
  else if (substr(pre,1,2)=="sp") pre <- "sphere"
  if (pre=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre=="sphere")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }
  
  Idr <- diag(d^r)
  RKr <- nu(r=r, diag(d))*2^(-d-r)*pi^(-d/2)
  if (!(pilot %in% c("dunconstr","dscalar")) | Sdr.flag)
  {
    Sd2r4 <- Sdr(d=d, r=2*r+4)
    if (nstage==2) Sd2r6 <- Sdr(d=d, r=2*r+6)
  }
  
  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  
  if (pilot=="unconstr")
  {
    ## psi4.mat is on data scale
    if (!Sdr.flag)
    {
       if (nstage==1)
        psi.fun <- (-1)^r*psifun1.unconstr(x=x, deriv.order=r, verbose=verbose, Sdr.flag=Sdr.flag, binned=FALSE)
       else if (nstage==2)
        psi.fun <- psifun2.unconstr(x=x, deriv.order=r, verbose=verbose, Sdr.flag=Sdr.flag, binned=FALSE, optim.fun=optim.fun)
    }
    else
    { 
      if (nstage==1)
        psi.fun <- (-1)^r*psifun1.unconstr(x=x, Sd2r4=Sd2r4, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
       else if (nstage==2)
        psi.fun <- psifun2.unconstr(x=x, Sd2r4=Sd2r4, Sd2r6=Sd2r6, binned=binned, bgridsize=bgridsize, deriv.order=r, verbose=verbose)
    }
    psi2r4.mat <- (-1)^r*invvec(psi.fun)   
    
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  else if (pilot=="dunconstr")
  {
    ## G2r4 is on data scale
    kfe.thin <- 1
    G2r4 <- Gdunconstr(x=x, d=d, r=r, n=n, nstage=nstage, verbose=verbose, Sdr.flag=Sdr.flag, optim.fun=optim.fun, thin=kfe.thin)
   
    vecPsi2r4 <- kfe(x=x, G=G2r4, binned=FALSE, deriv.order=2*r+4, deriv.vec=TRUE, add.index=FALSE, verbose=verbose, Sdr.flag=Sdr.flag, thin=kfe.thin)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  else if (pilot=="dscalar")
  {
    ## g2r4 is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
  }
  else
  {
    if (binned)
    {
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RKr)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      bin.par.star <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
    }
    
    ## psi4.mat is on pre-transformed data scale
    if (nstage==1)
      psi.fun <- psifun1(x.star, Sd2r4=Sd2r4, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    else if (nstage==2)
      psi.fun <- psifun2(x.star, Sd2r4=Sd2r4, Sd2r6=Sd2r6, pilot=pilot, binned=binned, bin.par=bin.par.star, deriv.order=r, verbose=verbose)$psir
    psi2r4.mat <- invvec(psi.fun)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  ## PI is estimate of AMISE
  pi.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    IdrvH <- Idr%x%vec(H)
    int.var <- 1/(det(H)^(1/2)*n)*nu(r=r, Hinv)*2^(-d-r)*pi^(-d/2)
    
    if (pilot=="dunconstr") 
    {
      ## eta functional form is not as accurate as kfe form
      ## eta2r4.hat <- eta.rs.kfe(x=x, r=4, s=2*r, A=H, B=diag(d), C=G2r4, verbose=verbose)
      ## pi.val <- int.var + (-1)^r*1/4*eta2r4.hat      
      pi.val <- drop(int.var + (-1)^r*1/4*vecPsi2r4 %*% (vec(diag(d^r) %x% vec(H) %x% vec(H))))  
    }
    else if (pilot=="dscalar")
    {
      eta2r4.hat <- eta.rs.kfe(x=x.star, r=4, s=2*r, A=H, B=diag(d), C=g2r4^2*diag(d), verbose=verbose)
      pi.val <- int.var + (-1)^r*1/4*eta2r4.hat
    }
    else
      pi.val <- int.var + (-1)^r*1/4* sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))

    pi.val <- drop(pi.val)
    return(pi.val)
  }

  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=vech(Hstart), f=pi.temp, print.level=2*as.numeric(verbose))    
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(vech(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  if (!(pilot %in% c("dunconstr","unconstr")))  H <- S12 %*% H %*% S12   ## back-transform
  
  if (!amise) return(H)
  else return(list(H = H, PI.star=amise.star))
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


Hpi.diag <- function(x, nstage=2, pilot="samse", pre="scale", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if(!is.matrix(x)) x <- as.matrix(x)
  if (substr(pre,1,2)=="sp") stop("Using pre-sphering won't give diagonal bandwidth matrix\n")
  if (substr(pilot,1,1)=="a") pilot <- "amse"
  else if (substr(pilot,1,1)=="s") pilot <- "samse"
  else if (substr(pilot,1,2)=="du") pilot <- "dunconstr"
  else if (substr(pilot,1,2)=="ds") pilot <- "dscalar"
 
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  RK <- (4*pi)^(-d/2)
  
  if (pilot=="amse" & (d>2 | r>0)) stop("SAMSE pilot selectors are better for higher dimensions and/or derivative r>0.")
  if (pilot=="samse" & r>0) stop("DSCALAR or DUNCONSTR pilot selectors are better for derivatives r>0.")
  if (pilot=="unconstr" | pilot=="dunconstr") stop("Unconstrained pilot selectors are not suitable for Hpi.diag.")
  
  if (substr(pre,1,2)=="sc") pre <- "scale"
  else if (substr(pre,1,2)=="sp") pre <- "sphere"
  if (pre=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre=="sphere")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }
  
  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par <- binning(x=x.star, bgridsize=bgridsize, H=diag(H.max)) 
  }
 
  Idr <- diag(d^r)
  ##RKr <- nu(r=r, diag(d))*2^(-d-r)*pi^(-d/2)
  
  if (pilot=="amse" | pilot=="samse")
  {
    ##Sd2r <- Sdr(d=d, r=2*r)
    Sd2r4 <- Sdr(d=d, r=2*r+4)
    if (nstage==2) Sd2r6 <- Sdr(d=d, r=2*r+6)
    
    if (nstage==1)
      psi.fun <- psifun1(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose, Sd2r4=Sd2r4)$psir
    else if (nstage==2)
      psi.fun <- psifun2(x.star, pilot=pilot, binned=binned, bin.par=bin.par, deriv.order=r, verbose=verbose, Sd2r4=Sd2r4, Sd2r6=Sd2r6)$psir
    psi2r4.mat <- invvec(psi.fun)
  }
  else if (pilot=="dscalar")
  {
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage)
  }
  
  if (d==2 & r==0 & (pilot=="amse" | pilot=="samse"))
  {
    ## diagonal bandwidth matrix for 2-dim has exact formula 
    psi40 <- psi.fun[1]
    psi22 <- psi.fun[6]
    psi04 <- psi.fun[16]
    s1 <- sd(x[,1])
    s2 <- sd(x[,2])
    h1 <- (psi04^(3/4)*RK/(psi40^(3/4)*(sqrt(psi40*psi04)+psi22)*n))^(1/6)
    h2 <- (psi40/psi04)^(1/4) * h1
    H <- diag(c(s1^2*h1^2, s2^2*h2^2))
    psimat4.D <- invvech(c(psi40, psi22, psi04))
    amise.star <- drop(n^(-1)*RK*(h1*h2)^(-1) + 1/4*c(h1,h2)^2 %*% psimat4.D %*% c(h1,h2)^2)
  }
  else
  {  
    ## PI is estimate of AMISE
    pi.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      Hinv <- chol2inv(chol(H))
      IdrvH <- Idr%x%vec(H)
      int.var <- 1/(det(H)^(1/2)*n)*nu(r=r, Hinv)*2^(-d-r)*pi^(-d/2)

      if (pilot=="dscalar")
      {
        eta2r4.hat <- eta.rs.kfe(x=x.star, r=4, s=2*r, A=H, B=diag(d), C=g2r4^2*diag(d), verbose=verbose)
        pi.val <- int.var + (-1)^r*1/4*eta2r4.hat
      }
      else
        pi.val <- int.var + (-1)^r*1/4* sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))
      
      return(drop(pi.val))
    }
    
    ## use normal reference bandwidth as initial condition
    if (missing(Hstart)) Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
    Hstart <- matrix.sqrt(Hstart)
    
    optim.fun1 <- tolower(substr(optim.fun,1,1))
    if (optim.fun1=="n")
    {
      result <- nlm(p=diag(Hstart), f=pi.temp, print.level=2*as.numeric(verbose))    
      H <- diag(result$estimate^2)
      amise.star <- result$minimum
    }
    else if (optim.fun1=="o")
    {  
      result <- optim(diag(Hstart), pi.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
      H <- diag(result$par) %*% diag(result$par)
      amise.star <- result$value
    }
    
    H <- S12 %*% H %*% S12  ## back-transform
  }

  if (!amise) return(H)
  else return(list(H = H, PI.star=amise.star))
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

lscv.mat <- function(x, H, binned=FALSE, bin.par, bgridsize, deriv.order=0, Sd2r, double.loop=FALSE, symm=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)
  
  if (!symm & !binned)
  {
    RK<-(4*pi)^(-d/2)
    ## fast version from J.E. Chacon 06/05/2011
    if (r==0)
    {
      Hinv<-chol2inv(chol(H))
      detH<-det(H)
      xH<-x%*%Hinv
      a<-rowSums(xH*x)
      M<-a%*%t(rep(1,n))+rep(1,n)%*%t(a)-2*(xH%*%t(x))
      m<-M[lower.tri(M)]
      em2<-exp(-m/2)
      lscv1<-(2*pi)^(-d/2)*2^(-d/2)*detH^(-1/2)*sum(sqrt(em2))
      lscv2<-(2*pi)^(-d/2)*detH^(-1/2)*sum(em2)        
      lscv <- RK/(n*detH^(1/2))+2*((1-1/n)*lscv1-2*lscv2)/(n*(n-1))
    }
    else if (r==1)
    {
      Hinv<-chol2inv(chol(H))
      H2inv<-Hinv%*%Hinv
      trHinv<-sum(diag(Hinv))
      detH<-det(H)
      
      xH<-x%*%Hinv
      a<-rowSums(xH*x)
      M<-a%*%t(rep(1,n))+rep(1,n)%*%t(a)-2*(xH%*%t(x))
      dv<-M[lower.tri(M)]
        
      xH2<-x%*%H2inv
      a2<-rowSums(xH2*x)
      M2<-a2%*%t(rep(1,n))+rep(1,n)%*%t(a2)-2*(xH2%*%t(x))
      dv2<-M2[lower.tri(M2)] 
      
      edv2<-exp(-dv/2)
      lscv1<-(2*pi)^(-d/2)*2^(-d/2)*detH^(-1/2)*sum(sqrt(edv2)*(dv2/4-trHinv/2))
      lscv2<-(2*pi)^(-d/2)*detH^(-1/2)*sum(edv2*(dv2-trHinv))
      lscv <- RK/(n*detH^(1/2))*trHinv*2^(-1)+(-1)*2*((1-1/n)*lscv1-2*lscv2)/(n*(n-1))
    }
    else if (r==2)
    {
      Hinv<-chol2inv(chol(H))
      H2inv<-Hinv%*%Hinv
      H3inv<-H2inv%*%Hinv
      trHinv<-sum(diag(Hinv))
      trH2inv<-sum(diag(H2inv))        
      detH<-det(H)
      
      xH<-x%*%Hinv
      a<-rowSums(xH*x)
      M<-a%*%t(rep(1,n))+rep(1,n)%*%t(a)-2*(xH%*%t(x))
      dv<-M[lower.tri(M)]
      
      xH2<-x%*%H2inv
      a2<-rowSums(xH2*x)
      M<-a2%*%t(rep(1,n))+rep(1,n)%*%t(a2)-2*(xH2%*%t(x))
      dv2<-M[lower.tri(M)]  
      
      xH3<-x%*%H3inv
      a3<-rowSums(xH3*x)
      M<-a3%*%t(rep(1,n))+rep(1,n)%*%t(a3)-2*(xH3%*%t(x))
      dv3<-M[lower.tri(M)]               
      
      edv2<-exp(-dv/2)
      lscv1<-(2*pi)^(-d/2)*2^(-d/2)*detH^(-1/2)*sum(sqrt(edv2)*(trH2inv/2-dv3/2+(-trHinv/2+dv2/4)^2))
      lscv2<-(2*pi)^(-d/2)*detH^(-1/2)*sum(edv2*(2*trH2inv-4*dv3+(-trHinv+dv2)^2))
      lscv <- RK/(n*detH^(1/2))*(trH2inv/2+trHinv^2/4)+2*((1-1/n)*lscv1-2*lscv2)/(n*(n-1))
    }
    else
    {
      ndifs <- n*(n-1)/2
      RK<-(4*pi)^(-d/2)        
      detH<-det(H)
      Hinv<-chol2inv(chol(H))
      
      xHinv<-x%*%Hinv
      xHinvx<-rowSums(xHinv*x)
      M<-xHinvx%*%t(rep(1,n))+rep(1,n)%*%t(xHinvx)-2*(xHinv%*%t(x)) #M=nxn matrix with m_ij=(X_i-X_j)^TH^{-1}(X_i-X_j)
      dv<-M[lower.tri(M)] # We only need the i<j values, I wonder if there is a way to obtain this not needing to calculate the full M
      
      edv2<-exp(-dv/2)
      P0<-Hinv
      
      kappas1<-matrix(nrow=ndifs,ncol=r)
      kappas2<-matrix(nrow=ndifs,ncol=r)
      kappas0<-numeric(r)        
      for (i in 1:r)
      {
        Hi1inv<-P0%*%Hinv
        trHi0inv<-sum(diag(P0))
        
        xHi1inv<-x%*%Hi1inv
        xHi1invx<-rowSums(xHi1inv*x)
        M<-xHi1invx%*%t(rep(1,n))+rep(1,n)%*%t(xHi1invx)-2*(xHi1inv%*%t(x))
        dvi1<-M[lower.tri(M)]
        
        kappas1[,i]<-(-2)^(i-1)*factorial(i-1)*(-trHi0inv/2^i+i*dvi1/2^(i+1))
        kappas2[,i]<-(-2)^(i-1)*factorial(i-1)*(-trHi0inv+i*dvi1)
        kappas0[i]<-(-2)^(i-1)*factorial(i-1)*(-trHi0inv/2^i)
        P0<-Hi1inv
      }
        
      nus1<-matrix(nrow=ndifs,ncol=r+1)
      nus2<-matrix(nrow=ndifs,ncol=r+1)        
      nus0<-numeric(r+1)         
      nus1[,1]<-1
      nus2[,1]<-1        
      nus0[1]<-1                
      for(j in 1:r)
      {
        js<-0:(j-1)
        if (j==1)
        {
          nus1[,2]<-kappas1[,1];nus2[,2]<-kappas2[,1];nus0[2]<-kappas0[1]
        }
        else
        {
          nus1[,j+1]<-rowSums(kappas1[,j:1]*nus1[,1:j]/matrix(rep(factorial(js)*factorial(rev(js)),ndifs),nrow=ndifs,byrow=TRUE))*factorial(j-1)
          nus2[,j+1]<-rowSums(kappas2[,j:1]*nus2[,1:j]/matrix(rep(factorial(js)*factorial(rev(js)),ndifs),nrow=ndifs,byrow=TRUE))*factorial(j-1)
          nus0[j+1]<-sum(kappas0[j:1]*nus0[1:j]/(factorial(js)*factorial(rev(js))))*factorial(j-1)
        }
      }
        
      lscv1<-(2*pi)^(-d/2)*2^(-d/2)*detH^(-1/2)*sum(sqrt(edv2)*nus1[,r+1])
      lscv2<-(2*pi)^(-d/2)*detH^(-1/2)*sum(edv2*nus2[,r+1])
      lscv<-(-1)^r*(RK/(n*detH^(1/2))*nus0[r+1]+2*((1-1/n)*lscv1-2*lscv2)/(n*(n-1)))
    }
  }
  else
  { 
    ## could implement these as recursive kfe 
    lscv1 <- kfe(x=x, G=2*H, inc=1, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, Sdr.mat=Sd2r, double.loop=double.loop, add.index=FALSE)
    lscv2 <- kfe(x=x, G=H, inc=0, binned=binned, bin.par=bin.par, bgridsize=bgridsize, deriv.order=2*r, Sdr.mat=Sd2r, double.loop=double.loop, add.index=FALSE)
    lscv <- drop(lscv1 - 2*lscv2)
    lscv <- (-1)^r*sum(vec(diag(d^r))*lscv)
  }

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

hlscv <- function(x, binned=TRUE, bgridsize, amise=FALSE, deriv.order=0)
{
  if (any(duplicated(x)))
    warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- length(x)
  d <- 1
  r <- deriv.order
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))

  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    bin.par <- binning(x, bgridsize=bgridsize, h=hnorm)
    lscv.1d.temp <- function(h) { return(lscv.1d(x=x, h=h, binned=binned, bin.par=bin.par, deriv.order=r)) }
  }
  else
  {
    if (r>0) stop("Unbinned hlscv not yet implemented for r>0.") 
    difs<-x%*%t(rep(1,n))-rep(1,n)%*%t(x)
    difs<-difs[lower.tri(difs)]  
    edifs<-exp(-difs^2/2)
    RK<-1/(2*sqrt(pi))

    lscv.1d.temp <- function(h)
    {
      lscv1 <- (1-1/n)*sum(edifs^(1/(2*h^2)))/(h*sqrt(2)*sqrt(2*pi))
      lscv2 <- 2*sum(edifs^(1/h^2))/(h*sqrt(2*pi))
      return(RK/(n*h)+2*(lscv1-lscv2)/(n^2-n))
    }    
  }
  opt <- optimise(f=lscv.1d.temp, interval=c(0.2*hnorm, 5*hnorm, tol=.Machine$double.eps))

  if (!amise) return(opt$minimum)
  else return(list(h=opt$minimum, LSCV=opt$objective))
}


Hlscv <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", trunc=4)
{
  if (any(duplicated(x))) warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order 
  ##RK <- (4*pi)^(-d/2)

  ## use normal reference selector as initial condn
  Hnorm <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  if (missing(Hstart)) Hstart <- matrix.sqrt(Hnorm)
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>4) binned <- FALSE
  if (binned) bin.par <- binning(x=x, H=diag(diag(Hnorm)))

  Sd2r <- Sdr(d=d, r=2*r)
  lscv.init <- lscv.mat(x=x, H=Hnorm, binned=FALSE, deriv.order=r, Sd2r=Sd2r, symm=FALSE)
  lscv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    lscv <- lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par, deriv.order=r, Sd2r=Sd2r, symm=FALSE)
    if (det(H) < 1/trunc*det(Hnorm) | det(H) > trunc*det(Hnorm) | abs(lscv) > trunc*abs(lscv.init)) lscv <- lscv.init 
    return(lscv)  
  }
    
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=vech(Hstart), f=lscv.mat.temp, print.level=2*as.numeric(verbose))   
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.opt <- result$minimum
  }
  else if (optim.fun1=="o")
  {
    result <- optim(vech(Hstart), lscv.mat.temp, control=list(trace=as.numeric(verbose)), method="Nelder-Mead")     
    H <- invvech(result$par) %*% invvech(result$par)
    amise.opt <- result$value
  }

  if (!amise) return(H)
  else return(list(H=H, LSCV=amise.opt))
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

Hlscv.diag <- function(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", trunc=4)
{
  if (any(duplicated(x))) warning("Data contain duplicated values: LSCV is not well-behaved in this case")
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  Hnorm <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  if (missing(Hstart)) Hstart <- matrix.sqrt(Hnorm)
  if (d>4) binned <- FALSE

  ## linear binning
  if (binned)
  {
    if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
    RK <- (4*pi)^(-d/2)
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    bin.par <- binning(x=x, bgridsize=bgridsize, H=sqrt(diag(diag(H.max))))
  }
  
  Sd2r <- Sdr(d=d, r=2*r)
  lscv.init <- lscv.mat(x=x, H=Hnorm, binned=FALSE, deriv.order=r, Sd2r=Sd2r, symm=FALSE)

  lscv.mat.temp <- function(diagH)
  {
    H <- diag(diagH^2)
    lscv <- lscv.mat(x=x, H=H, binned=binned, bin.par=bin.par, deriv.order=r, Sd2r=Sd2r, symm=FALSE)
    if (det(H) < 1/trunc*det(Hnorm) | det(H) > trunc*det(Hnorm) | abs(lscv) > trunc*abs(lscv.init)) lscv <- lscv.init 
    return(lscv)  
  }
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=diag(Hstart), f=lscv.mat.temp, print.level=2*as.numeric(verbose))
    H <- diag(result$estimate^2)
    amise.opt <- result$minimum
  }
  else if (optim.fun1=="o")
  {
    result <- optim(diag(Hstart), lscv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par^2)
    amise.opt <- result$value
  }
  
  if (!amise) return(H)
  else return(list(H=H, LSCV=amise.opt))
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

Hbcv <- function(x, whichbcv=1, Hstart, amise=FALSE, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)

  ## use normal reference b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- Hmax
  if (missing(Hstart)) Hstart <- matrix.sqrt(0.9*Hmax)

  bcv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(bcv.mat(x, H, whichbcv*H)$bcv)
  } 
 
  result <- optim(vech(Hstart), bcv.mat.temp, method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)), lower=-vech(matrix.sqrt(up.bound)), control=list(trace=as.numeric(verbose)))
  H <- invvech(result$par) %*% invvech(result$par)
  amise.opt <- result$value
  
  if (!amise) return(H)
  else return(list(H = H, BCV=amise.opt))
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

Hbcv.diag <- function(x, whichbcv=1, Hstart, amise=FALSE, verbose=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  ## use maximally smoothed b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- diag(Hmax)
  
  if (missing(Hstart)) Hstart <- 0.9*matrix.sqrt(Hmax)

  bcv.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH) 
    return(bcv.mat(x, H, whichbcv*H)$bcv)
  } 
 
  result <- optim(diag(Hstart), bcv.mat.temp, method="L-BFGS-B", upper=sqrt(up.bound), control=list(trace=as.numeric(verbose)))
  H <- diag(result$par) %*% diag(result$par)
  amise.opt <- result$value
  if (!amise) return(H)
  else return(list(H = H, BCV=amise.opt))
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

gamse.scv <- function(x.star, Sd6, d, Sigma.star, Hamise, n, binned=FALSE, bin.par, bgridsize, verbose=FALSE, nstage=1, Theta6=FALSE)
{
  if (nstage==0)
  {
    psihat6.star <- psins(r=6, Sigma=Sigma.star, deriv.vec=TRUE, Sdr.mat=Sd6) 
  }
  else if (nstage==1)
  {  
    g6.star <- gsamse(Sigma.star, n, 6) 
    G6.star <- g6.star^2 * diag(d)
    if (Theta6) psihat6.star <- kfe(x=x.star, bin.par=bin.par, Sdr.mat=Sd6, deriv.order=6, G=G6.star, deriv.vec=FALSE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    else psihat6.star <- kfe(x=x.star, bin.par=bin.par, Sdr.mat=Sd6, deriv.order=6, G=G6.star, deriv.vec=TRUE, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }

  if (Theta6)
  {
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
    gamse <- (num/(den*n))^(1/(d+6)) 
  }
  else
  {  
    ## updated constants using Chacon & Duong (2010) notation
    Cmu1Cmu1 <- drop(1/4*psihat6.star %*% (Hamise %x% diag(d^4) %x% Hamise) %*% psihat6.star)
    Cmu1Cmu2 <- 3/4*(4*pi)^(-d/2)*drop(vec(Hamise %x% diag(d) %x% Hamise) %*% psihat6.star)
    Cmu2Cmu2 <- 1/64*(4*pi)^(-d)*(4*tr(Hamise%*%Hamise) + (d+8)*tr(Hamise)^2)
    num <- 2 * (d+4) * Cmu2Cmu2
    den <- -(d+2) * Cmu1Cmu2 + sqrt((d+2)^2 * Cmu1Cmu2^2 + 8*(d+4)*Cmu1Cmu1 * Cmu2Cmu2)
    gamse <- (num/(den*n))^(1/(d+6)) 
  }
  return(gamse)
}


###############################################################################
# Estimate unconstrained G_AMSE pilot bandwidth for SCV for 2 to 6 dim
# (J.E. Chacon)
#
# Returns
# G_AMSE pilot bandwidth
###############################################################################

Gunconstr.scv <- function(x, Sd6, Sd4, binned=FALSE, bin.par, bgridsize, rel.tol=10^-10, verbose=FALSE, nstage=1, Sdr.flag=TRUE, optim.fun="nlm")
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  ## stage 1 of plug-in
  if (nstage==1)
  {  
    G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
    psihat6 <- kfe(x=x, deriv.order=6, Sdr.mat=Sd6, G=G6, deriv.vec=TRUE, add.index=FALSE, binned=binned, bin.par=bin.par, bgridsize=bgridsize, verbose=verbose, Sdr.flag=Sdr.flag)
  }
  else if (nstage==0)
  {
    psihat6 <- psins(r=6, Sigma=S, deriv.vec=TRUE, Sdr.mat=Sd6, Sdr.flag=Sdr.flag) 
  }
  
  ## constants for normal reference
  ##D4phi0 <- drop(dmvnorm.deriv(x=rep(0,d), deriv.order=4, Sdr.mat=Sd4))
  D4phi0 <- DrL0(d=d,r=4, Sdr.mat=Sd4, Sdr.flag=Sdr.flag)
  Id4 <- diag(d^4)

  ## asymptotic squared bias for r = 4
  AB2<-function(vechG){
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(Kpow(A=Ginv12,pow=4)%*%D4phi0)*2^(-(d+4)/2) + (t(vec(G))%x%Id4)%*%psihat6
    return (sum(AB^2))
  }

  Hstart <- (4/(d+2))^(2/(d+4))*n^(-2/(d+4))*S
  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
     result <- nlm(p=vech(Hstart), f=AB2, print.level=2*as.logical(verbose))    
     G4 <- result$estimate
  }
  else    
  { 
    result <- optim(vech(matrix.sqrt(Hstart)), AB2, method="BFGS", control=list(trace=as.numeric(verbose)))
    G4 <- result$par
  }   
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
  scv1 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(2*h^2+2*g^2), binned=binned, inc=inc)
  scv2 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(h^2+2*g^2), binned=binned, inc=inc)
  scv3 <- kfe.1d(x=x, deriv.order=2*r, bin.par=bin.par, g=sqrt(2*g^2), binned=binned, inc=inc)

  bias2 <- (-1)^r*(scv1 - 2*scv2 + scv3)
  if (bias2 < 0) bias2 <- 0
  scv <- (n*h)^(-1)*(4*pi)^(-1/2)*2^(-r)*OF(2*r) + bias2

  return(scv)
}

scv.mat <- function(x, H, G, binned=FALSE, bin.par, bgridsize, verbose=FALSE, deriv.order=0, Sd2r, symm=FALSE, Sdr.flag=FALSE)
{
  d <- ncol(x)
  n <- nrow(x)
  r <- deriv.order
  vId <- vec(diag(d))
  Hinv <- chol2inv(chol(H))

  if (!symm)
  {
    scv1 <- eta.kfe.y(x=x, G=2*H+2*G, deriv.order=2*r)
    scv2 <- eta.kfe.y(x=x, G=H+2*G, deriv.order=2*r)
    scv3 <- eta.kfe.y(x=x, G=2*G, deriv.order=2*r)
    bias2 <- (-1)^r*(scv1 - 2*scv2 + scv3)
    if (bias2 < 0) bias2 <- 0
  }
  else
  {
    scv1 <- kfe(x=x, G=2*H + 2*G, deriv.order=2*r, Sdr.mat=Sd2r, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)$psir
    scv2 <- kfe(x=x, G=H + 2*G, deriv.order=2*r, Sdr.mat=Sd2r, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)$psir
    scv3 <- kfe(x=x, G=2*G, deriv.order=2*r, Sdr.mat=Sd2r, inc=1, bin.par=bin.par, binned=binned, bgridsize=bgridsize, verbose=verbose, add.index=FALSE)$psir
    
    bias2 <- drop((-1)^r*Kpow(vId,r) %*% (scv1 - 2*scv2 + scv3))
    if (bias2 < 0) bias2 <- 0
  }
  scvmat <- 1/(det(H)^(1/2)*n)*nu(r=r, Hinv)*2^(-d-r)*pi^(-d/2) + bias2
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


Hscv <- function(x, nstage=2, pre="sphere", pilot="samse", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm", Sdr.flag=TRUE)
{
  n <- nrow(x)
  d <- ncol(x)
  ##S <- var(x)
  r <- deriv.order

  if(!is.matrix(x)) x <- as.matrix(x)
  if (substr(pilot,1,1)=="a") pilot <- "amse"
  else if (substr(pilot,1,1)=="s") pilot <- "samse"
  else if (substr(pilot,1,1)=="u") pilot <- "unconstr"
  else if (substr(pilot,1,2)=="du") pilot <- "dunconstr"                     
  else if (substr(pilot,1,2)=="ds") pilot <- "dscalar"

  if (pilot=="amse" & (d>2 | r>0)) stop("AMSE pilot selectors not defined for d>2 and/or r>0.")
  if ((pilot=="samse" | pilot=="unconstr") & r>0) stop("DSCALAR or DUNCONSTR pilot selectors are better for derivatives r>0.")
    
  if (substr(pre,1,2)=="sc") pre <- "scale"
  else if (substr(pre,1,2)=="sp") pre <- "sphere"
  
  if (pre=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre=="sphere")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }
    
  Sd2r <- Sdr(d=d,r=2*r)
  ##vId <- vec(diag(d))
  ##RKr <- nu(r=r, diag(d))*2^(-d-r)*pi^(-d/2)
  RK <- (4*pi)^(-d/2)

  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  
  if (binned)
  {
    if (pilot=="unconstr" | pilot=="dunconstr")
    {
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
      bin.par <- binning(x=x, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
    }
    else
    {  
      H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      bin.par <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
    }
  }
  if (pilot=="unconstr")
  {
    ## Gu pilot matrix is on data scale
    Sd4 <- Sdr(d=d, r=4)
    Sd6 <- Sdr(d=d, r=6)
    Gu <- Gunconstr.scv(x=x, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, Sd4=Sd4, nstage=nstage-1, Sdr.flag=Sdr.flag, optim.fun=optim.fun)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  else if (pilot=="dunconstr")
  {
    ## Gu pilot matrix is on data scale
    Gu <- Gdunconstr(x=x, d=d, r=r, n=n, nstage=nstage, verbose=verbose, scv=TRUE, Sdr.flag=Sdr.flag, optim.fun=optim.fun)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x)
  }
  else if (pilot=="dscalar")
  {
    ## Gs is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, d=d, r=r, n=n, nstage=nstage, verbose=verbose, scv=TRUE)
    Gs <- g2r4^2*diag(d)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
  }
  else
  {
    ## Gs is on transformed data scale
    
    Hamise <- Hpi(x=x.star, nstage=1, deriv.order=r, pilot=pilot, pre="sphere", binned=TRUE, bgridsize=bgridsize, verbose=verbose, Sdr.flag=Sdr.flag, optim.fun=optim.fun) 
    if (any(is.na(Hamise)))
    {
      warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
      Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    }
    Sd4 <- Sdr(d=d, r=4)
    Sd6 <- Sdr(d=d, r=6)
    gs <- gamse.scv(x.star=x.star, d=d, Sigma.star=var(x.star), Hamise=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, nstage=nstage-1)
    Gs <- gs^2*diag(d)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  ## SCV is estimate of AMISE
  scv.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    if (pilot=="samse" | pilot=="amse" | pilot=="dscalar"){ Gpilot <- Gs; xx <- x.star }
    else if (pilot=="unconstr" | pilot=="dunconstr") { Gpilot <- Gu; xx <- x }
    return(scv.mat(x=xx, H=H, G=Gpilot, binned=binned, bin.par=bin.par, verbose=FALSE, Sd2r=Sd2r, deriv.order=r, symm=FALSE))
  }

  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=vech(Hstart), f=scv.mat.temp, print.level=2*as.numeric(verbose))    
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(vech(Hstart), scv.mat.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  if (!(pilot %in% c("dunconstr","unconstr")))  H <- S12 %*% H %*% S12   ## back-transform

  if (!amise) return(H)
  else return(list(H = H, SCV.star=amise.star))
}


Hscv.diag <- function(x, nstage=2, pre="scale", pilot="samse", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  RK <- (4*pi)^(-d/2)
  
  if(!is.matrix(x)) x <- as.matrix(x)
  if (substr(pilot,1,1)=="a") pilot <- "amse"
  else if (substr(pilot,1,1)=="s") pilot <- "samse"
  else if (substr(pilot,1,2)=="du") pilot <- "dunconstr"
  else if (substr(pilot,1,2)=="ds") pilot <- "dscalar"

  if (pilot=="amse" & (d>2 | r>0)) stop("SAMSE pilot selectors are better for higher dimensions and/or derivative r>0.")
  if (pilot=="samse" & r>0) stop("DSCALAR pilot selectors are better for derivatives r>0.")
  if (pilot=="unconstr" | pilot=="dunconstr") stop("Unconstrained pilot selectors are not suitable for Hscv.diag.")

  if (substr(pre,1,2)=="sc") pre <- "scale"
  else if (substr(pre,1,2)=="sp") pre <- "sphere"
  if (pre=="sphere") stop("Using pre-sphering doesn't give a diagonal bandwidth matrix\n")

  if (pre=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre=="sphere")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }

  Sd2r <- Sdr(d=d,r=2*r)
  if (d > 4) binned <- FALSE
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  if (d>=4 & nstage==2) bgridsize <- rep(11,d)
  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    bin.par.star <- binning(x=x.star, bgridsize=bgridsize, H=sqrt(diag(H.max))) 
  }
  
  if (pilot=="dscalar")
  {
    ## Gs is on pre-transformed data scale
    g2r4 <- gdscalar(x=x.star, r=r, n=n, d=d, verbose=verbose, nstage=nstage, scv=TRUE)
    Gs <- g2r4^2*diag(d)
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
  }
  else
  {
    ## Gs is on transformed data scale
      Hamise <- Hpi(x=x.star, nstage=1, pilot=pilot, pre="sphere", binned=binned, bgridsize=bgridsize, verbose=verbose, optim.fun=optim.fun) 
      if (any(is.na(Hamise)))
      {
        warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
        Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
      }
      ##Sd4 <- Sdr(d=d, r=4)
      Sd6 <- Sdr(d=d, r=6)
      gs <- gamse.scv(x.star=x.star, d=d, Sigma.star=var(x.star), Hamise=Hamise, n=n, binned=binned, bgridsize=bgridsize, verbose=verbose, Sd6=Sd6, nstage=nstage-1)
    
    Gs <- gs^2*diag(d)

    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) Hstart <- (4/(n*(d+2*r+2)))^(2/(d+2*r+4)) * var(x.star)
    else Hstart <- Sinv12 %*% Hstart %*% Sinv12
  }

  scv.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(scv.mat(x.star, H, Gs, binned=binned, bin.par=bin.par.star, verbose=FALSE, Sd2r=Sd2r, deriv.order=r, symm=FALSE))
  }

  Hstart <- matrix.sqrt(Hstart)

  ## back-transform
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=diag(Hstart), f=scv.mat.temp, print.level=2*as.numeric(verbose))    
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="o")
  {  
    result <- optim(diag(Hstart), scv.mat.temp, method="Nelder-Mead", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }
  H <- S12 %*% H %*% S12

  if (!amise) return(H)
  else return(list(H = H, SCV.star=amise.star))
}



