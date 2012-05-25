###############################################################################
## Estimate plug-in unconstrained bandwidth for psi_0 (for 2-sample test)
##
## Returns
## Plug-in bandwidth
###############################################################################


hpi.kfe <- function(x, nstage=2, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0)
{
  n <- length(x)
  d <- 1
  r <- deriv.order
  k <- 2 ## kernel order
  Kr0 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=r)    
  mu2K <- 1  
  
  if (nstage==2)
  {
    psi4.hat <- psins.1d(r=r+k+2, sigma=sd(x))
    gamse2 <- (factorial(r+2)*Kr0/(mu2K*psi4.hat*n))^(1/(r+k+3))
    psi2.hat <- kfe.1d(x=x, g=gamse2, deriv.order=r+k, inc=1, binned=binned)
  }
  else 
  {
    psi2.hat <- psins.1d(r=r+k, sigma=sd(x))
  }

  ## formula for bias annihliating bandwidths from Wand & Jones (1995, p.70)
  gamse <- (factorial(r)*Kr0/(-mu2K*psi2.hat*n))^(1/(r+k+1)) 
  
  return(gamse)
}


Hpi.kfe <- function(x, nstage=2, pilot="dscalar", pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("Currently only deriv.order=0 is implemented")
   
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order

  if(!is.matrix(x)) x <- as.matrix(x)
  if (substr(pilot,1,2)=="du") pilot <- "dunconstr"                     
  else if (substr(pilot,1,2)=="ds") pilot <- "dscalar"
 
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

  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)

  if (nstage==2)
  {  
    ## stage 1
    if (pilot=="dscalar")
    {
      psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
      ##h2 <- gdscalar(x=x.star, d=d, r=r-2, n=n, verbose=verbose, nstage=nstage, scv=FALSE)     
      ## gdscalar not used because it's an implicit computation without
      ## symmetriser matrices and is slower than direct computation with 
      ## symmetriser matrices. 

      A1 <- drop(t(D2K0) %*% D2K0)
      A2 <- drop(t(D2K0) %*% t(vec(diag(d)) %x% diag(d^2)) %*% psi4.ns) 
      A3 <- drop(t(psi4.ns) %*% (vec(diag(d)) %*% t(vec(diag(d))) %x% diag(d^2)) %*% psi4.ns)

      ## Special case from Chacon & Duong (2010): bias minimisation
      ##h2 <- ((4*d+8)*A1/(-d*A2 + sqrt(d^2*A2^2 + (8*d+16)*A1*A3)))^(1/(d+4))*n^(-1/(d+4))
      ## Special from Chacon & Duong (2010): bias annihilation
      h2 <- (-A1/(2*A2*n))^(1/(d+4))
      H2 <- h2^2*diag(d)   
      psi2.hat <- kfe(x=x.star, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
    else if (pilot=="dunconstr")
    {
      psi4.ns <- psins(r=r+4, Sigma=var(x), deriv.vec=TRUE)
 
      amse2.temp <- function(vechH)
      { 
        H <- invvech(vechH) %*% invvech(vechH)
        Hinv <- chol2inv(chol(H))
        Hinv12 <- matrix.sqrt(Hinv)
        amse2.temp <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
        return(sum((amse2.temp)^2)) 
      }
      
      Hstart2 <- matrix.sqrt(Gns(r=r+2, n=n, Sigma=var(x)))
      optim.fun1 <- tolower(substr(optim.fun,1,1))
 
      if (optim.fun1=="n")
      {
         result <- nlm(p=vech(Hstart2), f=amse2.temp, print.level=2*as.numeric(verbose))    
         H2 <- invvech(result$estimate) %*% invvech(result$estimate)
      }
      else
      {
         result <- optim(vech(Hstart2), amse2.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
         H2 <- invvech(result$par) %*% invvech(result$par)
      }
 
      psi2.hat <- kfe(x=x, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
  }
  else
  {
    if (pilot=="dscalar") psi2.hat <- psins(r=r+2, Sigma=var(x.star), deriv.vec=TRUE)
    else if (pilot=="dunconstr") psi2.hat <- psins(r=r+2, Sigma=var(x), deriv.vec=TRUE)    
  }

  if (pilot=="dscalar") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x.star))}
  else if (pilot=="dunconstr") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x))}
  
  ## stage 2
  amse.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    amse.val <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.val^2))) 
  }
  
  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=vech(Hstart), f=amse.temp, print.level=2*as.numeric(verbose)) 
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(vech(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  
  ## back-transform
  if (pilot=="dscalar") H <- S12 %*% H %*% S12 

  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}


Hpi.diag.kfe <- function(x, nstage=2, pilot="dscalar", pre="scale", Hstart, deriv.order=0, binned=FALSE, bgridsize, amise=FALSE, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("Currently only deriv.order=0 is implemented")

  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)

  if (substr(pre,1,2)=="sp") stop("Using pre-sphering won't give diagonal bandwidth matrix\n")
  if (substr(pilot,1,2)=="ds") pilot <- "dscalar"
  
  if (substr(pre,1,2)=="sc") pre <- "scale"
  if (pre=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }

  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
   
    if (pilot=="dscalar")
    {
      ## h2 is on pre-transformed data scale
      A1 <- drop(t(D2K0) %*% D2K0)
      A2 <- drop(t(D2K0) %*% t(vec(diag(d)) %x% diag(d^2)) %*% psi4.ns) 
      A3 <- drop(t(psi4.ns) %*% (vec(diag(d)) %*% t(vec(diag(d))) %x% diag(d^2)) %*% psi4.ns)
      h2 <- ((4*d+8)*A1/(-d*A2 + sqrt(d^2*A2^2 + (8*d+16)*A1*A3)))^(1/(d+4))*n^(-1/(d+4))
      H2 <- h2^2*diag(d)
    }  
    psi2.hat <- kfe(x=x.star, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose) 
  }
  else
    psi2.hat <- psins(r=r+2, Sigma=var(x.star), deriv.vec=TRUE)
 
  ## stage 2
  amse.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    amse.val <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.val^2))) 
  }

  if (missing(Hstart)) Hstart <- Hns(x=x.star, deriv.order=r)
  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- tolower(substr(optim.fun,1,1))
  if (optim.fun1=="n")
  {
    result <- nlm(p=diag(Hstart), f=amse.temp, print.level=2*as.numeric(verbose))    
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(diag(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }
  ## back-transform
  if (pilot=="dscalar") H <- S12 %*% H %*% S12 
  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}


##############################################################################
## Test statistic for multivariate 2-sample test
##############################################################################

kde.test <- function(x1, x2, H1, H2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, binned=FALSE, bgridsize, verbose=FALSE, pilot="dscalar")
{
  ##if (binned) warning("From ks 1.8.8, binned estimation now only applies to the calculation of the bandwidths H1 and H2, and not the pvalue.") 
  
  if (is.vector(x1) & is.vector(x2)) 
    return(kde.test.1d(x1=x1, x2=x2, h1=h1, h2=h2, psi1=psi1, psi2=psi2, var.fhat1=var.fhat1, var.fhat2=var.fhat2, binned=binned, bgridsize=bgridsize, verbose=verbose))
    
  if (!is.matrix(x1)) x1 <- as.matrix(x1)
  if (!is.matrix(x2)) x2 <- as.matrix(x2)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  d <- ncol(x1)
  K0 <- drop(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0))
  
  ## kernel estimation for components of test statistic
  if (missing(H1)) H1 <- Hpi.kfe(x1, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize, verbose=verbose, pilot=pilot)
  if (missing(H2)) H2 <- Hpi.kfe(x2, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize, verbose=verbose, pilot=pilot)

  symm <- FALSE ## don't use symmetriser matrices in psi functional calculations
  if (missing(psi1)) psi1 <- eta.kfe.y(x=x1, y=x1, G=H1, verbose=verbose, symm=symm)      
  if (missing(psi2)) psi2 <- eta.kfe.y(x=x2, y=x2, G=H2, verbose=verbose, symm=symm)
 
  S1 <- var(x1)
  S2 <- var(x2)
    
  if (missing(var.fhat1))
  {
    H1.r1 <- Hns(x=x1, deriv.order=1)
    fhat1.r1 <- kdde(x=x1, H=H1.r1, deriv.order=1, eval.points=apply(x1, 2, mean))$estimate
    var.fhat1 <- drop(fhat1.r1 %*% S1 %*% t(fhat1.r1))
  }
  psi12 <- eta.kfe.y(x=x1, G=H1, y=x2, verbose=verbose, symm=symm) 
  
  if (missing(var.fhat2))
  {
    H2.r1 <- Hns(x=x2, deriv.order=1)
    fhat2.r1 <- kdde(x=x2, H=H2.r1, deriv.order=1, eval.points=apply(x2, 2, mean))$estimate
    var.fhat2 <- drop(fhat2.r1 %*% S2 %*% t(fhat2.r1))
  }
  psi21 <- eta.kfe.y(x=x2, G=H2, y=x1, verbose=verbose, symm=symm)
  
  ## test statistic + its parameters
  
  T.hat <- drop(psi1 + psi2 - (psi12 + psi21))
  muT.hat <- (n1^(-1)*det(H1)^(-1/2) + n2^(-1)*det(H2)^(-1/2))*K0
  varT.hat <- 3*(n1*var.fhat1 + n2*var.fhat2)/(n1+n2) *(1/n1+1/n2) 

  zstat <- (T.hat-muT.hat)/sqrt(varT.hat)
  pval <- 1-pnorm(zstat)
  if (pval==0) pval <- pnorm(-abs(zstat)) 
 
  val <- list(Tstat=T.hat, zstat=zstat, pvalue=pval, mean=muT.hat, var=varT.hat, var.fhat1=var.fhat1, var.fhat2=var.fhat2, n1=n1, n2=n2, H1=H1, H2=H2, psi1=psi1, psi12=psi12, psi21=psi21, psi2=psi2)
  return(val)
}     


kde.test.1d <- function(x1, x2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, binned=FALSE, bgridsize, verbose=FALSE)
{
  n1 <- length(x1)
  n2 <- length(x2)
  d <- 1
  K0 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=0)

  s1 <- sd(x1)
  s2 <- sd(x2)

  ## kernel estimation for components of test statistic
  if (missing(h1)) h1 <- hpi.kfe(x1, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize)
  if (missing(h2)) h2 <- hpi.kfe(x2, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize)

  if (missing(psi1)) psi1 <- eta.kfe.y.1d(x=x1, y=x1, g=h1, verbose=verbose)
  if (missing(psi2)) psi2 <- eta.kfe.y.1d(x=x2, y=x2, g=h2, verbose=verbose)
   
  if (missing(var.fhat1))
  {
    h1.r1 <- hns(x=x1, deriv.order=1) 
    fhat1.r1 <- kdde(x=x1, h=h1.r1, deriv.order=1, eval.points=mean(x1))$estimate
    var.fhat1 <- fhat1.r1^2*s1^2
  }
  psi12 <- eta.kfe.y.1d(x=x1, g=h1, y=x2, verbose=verbose)
  
  if (missing(var.fhat2))
  {
    h2.r1 <- hns(x=x2, deriv.order=1) 
    fhat2.r1 <- kdde(x=x2, h=h2.r1, deriv.order=1, eval.points=mean(x2))$estimate
    var.fhat2 <- fhat2.r1^2*s2^2
  }
  psi21 <- eta.kfe.y.1d(x=x2, g=h2, y=x1, verbose=verbose)
 

  ## test statistic + its parameters
  T.hat <- drop(psi1 + psi2 - (psi12 + psi21))
  muT.hat <- ((n1*h1)^(-1) + (n2*h2)^(-1))*K0
  varT.hat <- 3*(n1*var.fhat1 + n2*var.fhat2)/(n1+n2) *(1/n1+1/n2) 
  zstat <- (T.hat-muT.hat)/sqrt(varT.hat)
  pval <- 1-pnorm(zstat)
  if (pval==0) pval <- pnorm(-abs(zstat)) 
 
  val <- list(Tstat=T.hat, zstat=zstat, pvalue=pval, mean=muT.hat, var=varT.hat, var.fhat1=var.fhat1, var.fhat2=var.fhat2, n1=n1, n2=n2, h1=h1, h2=h2, psi1=psi1, psi12=psi12, psi21=psi21, psi2=psi2)
  return(val)

}


