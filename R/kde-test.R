###############################################################################
## Estimate plug-in unconstrained bandwidth for psi_0 (for 2-sample test)
##
## Returns
## Plug-in bandwidth
###############################################################################

Hpi.kfe <- function(x, nstage=2, Hstart, deriv.order=0, double.loop=FALSE, amise=FALSE, verbose=FALSE)
{
  if (deriv.order!=0) stop("Currently only deriv.order=0 is implemented")
  
  n <- nrow(x)
  d <- ncol(x)

  ## use normal reference bandwidth as initial condition 
  if (missing(Hstart)) { r <- 4; Hstart <- 2*(2/(n*(d+r)))^(2/(d+r+2)) * var(x) }
  Hstart <- matrix.sqrt(Hstart)
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)
  m2K2 <- 1/2*(4*pi)^(-d/2)*vec(diag(d))
  ##m0K2 <- (4*pi)^(-d/2)
  
  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=4, Sigma=var(x), deriv.vec=TRUE)
    amse2.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      Hinv <- chol2inv(chol(H))
      Hinv12 <- matrix.sqrt(Hinv)
      amse2.temp <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
      return(sum((amse2.temp)^2)) 
    }
    result <- optim(vech(Hstart), amse2.temp, method="BFGS")
    H2 <- invvech(result$par) %*% invvech(result$par)
    
    psi2.hat <- kfe(x=x, G=H2, deriv.order=2, double.loop=double.loop, add.index=FALSE, verbose=verbose)
  }
  else
    psi2.hat <- psins(r=2, Sigma=var(x), deriv.vec=TRUE)

  ## stage 2
  amse.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    amse.temp <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.temp^2))) 
  }
  r <- 2; Hstart <- 2*(2/(n*(d+r)))^(2/(d+r+2)) * var(x)
  Hstart <- matrix.sqrt(Hstart)
  result <- optim(vech(Hstart), amse.temp, method="BFGS")
  H <- invvech(result$par) %*% invvech(result$par)

  if (!amise)
    return(H)
  else
    return(list(H = H, PI=result$value))
}


Hpi.diag.kfe <- function(x, nstage=2, Hstart, deriv.order=0, binned=FALSE, double.loop=FALSE, amise=FALSE, verbose=FALSE)
{
  if (deriv.order!=0) stop("Currently only deriv.order=0 is implemented")

  n <- nrow(x)
  d <- ncol(x)

  ## use normal reference bandwidth as initial condition 
  if (missing(Hstart)) { r <- 4; Hstart <- 2*(2/(n*(d+r)))^(2/(d+r+2)) * var(x) }
  Hstart <- matrix.sqrt(Hstart)
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)
  m2K2 <- 1/2*(4*pi)^(-d/2)*vec(diag(d))

  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=4, Sigma=var(x), deriv.vec=TRUE)
    
    amse2.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      Hinv <- chol2inv(chol(H))
      Hinv12 <- matrix.sqrt(Hinv)
      amse2.temp <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
      return(sum((amse2.temp)^2)) 
    }
    result <- optim(diag(Hstart), amse2.temp, method="BFGS")
    H2 <- diag(result$par) %*% diag(result$par)

    
    if (binned) bin.par <- binning(x, H=H2)
    psi2.hat <- kfe(x=x, G=H2, deriv.order=2, double.loop=double.loop, add.index=FALSE, binned=binned, bin.par=bin.par, verbose=verbose) 
  }
  else
    psi2.hat <- psins(r=2, Sigma=var(x), deriv.vec=TRUE)
 
  ## stage 2
  amse.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    amse.temp <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.temp^2))) 
  }
  r <- 2; Hstart <- 2*(2/(n*(d+r)))^(2/(d+r+2)) * var(x)
  Hstart <- matrix.sqrt(Hstart)
  result <- optim(diag(Hstart), amse.temp, method="BFGS")
  H <- diag(result$par) %*% diag(result$par)
  
  if (!amise)
    return(H)
  else
    return(list(H = H, PI=result$value))
}


#######################################################################################################
## Kernel estimator of xi = int f(x)^3 dx = int f(x)^2 f(x) dx
#######################################################################################################

xi <- function(x, G)
{
  n <- nrow(x)
  d <- ncol(x)
  difs <- differences(x, upper=FALSE)
  K12 <- dmvnorm(difs, mean=rep(0,d), sigma=G)
  K12K13 <- 0

  for (i in 1:n)
  {
    K12K13 <- K12K13 + sum(K12[((i-1)*n+1):(i*n)] %x% K12[((i-1)*n+1):(i*n)])
    ##cat(i)
  }

  return(K12K13/n^3)
}




#######################################################################################################
## Test statistic for 2-sample test
#######################################################################################################

kde.test <- function(x1, x2, H1, H2, psi1, psi2, fhat1, fhat2, var.fhat1, var.fhat2, double.loop=FALSE, binned=FALSE, verbose=FALSE)
{
  n1 <- nrow(x1)
  n2 <- nrow(x2)

  d <- ncol(x1)
  n <- nrow(x1)
  K0 <- drop(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0))
 
  ## kernel estimation for components of test statistic
  if (missing(H1))
    if (binned) H1 <- Hpi.diag.kfe(x1, nstage=2, double.loop=double.loop, deriv.order=0, binned=TRUE, verbose=verbose)
    else H1 <- Hpi.kfe(x1, nstage=2, double.loop=double.loop, deriv.order=0, verbose=verbose)
  if (missing(H2))
    if (binned) H2 <- Hpi.diag.kfe(x2, nstage=2, double.loop=double.loop, deriv.order=0, binned=TRUE, verbose=verbose)
    else H2 <- Hpi.kfe(x2, nstage=2, double.loop=double.loop, deriv.order=0, verbose=verbose)

  if (missing(psi1)) psi1 <- kfe(x=x1, G=H1, deriv.order=0, double.loop=double.loop, add.index=FALSE, binned=binned, verbose=verbose)
  if (missing(psi2)) psi2 <- kfe(x=x2, G=H2, deriv.order=0, double.loop=double.loop, add.index=FALSE, binned=binned, verbose=verbose)
  
  if (!missing(fhat1))
  {
    fhat1 <- find.nearest.gridpts(x=rbind(x1,x2), gridx=fhat1$eval.points, f=fhat1$estimate)$fx
    psi12 <- sum(tail(fhat1, n=n2))/n2
    var.fhat1 <- var(head(fhat1, n=n1))
  }
  else
  {  
    if (missing(var.fhat1))
    {
      #fhat1.x1.params <- kde.points.sum(x=x1, H=H1, eval.points=x1)
      #var.fhat1 <- (fhat1.x1.params$sumsq - fhat1.x1.params$sum^2/n1)/(n1-1)
      S1 <- var(x1)
      H1.r1 <- Hamise.mixt(mus=rep(0,d), Sigmas=S1, samp=n1, props=1, deriv.order=1)
      fhat1.r1 <- kdde(x=x1, H=H1.r1, deriv.order=1, eval.points=apply(x1, 2, mean))$estimate
      var.fhat1 <- drop(fhat1.r1 %*% S1 %*% t(fhat1.r1))
    }
    psi12 <- kde.points.sum(x=x1, H=H1, eval.points=x2, verbose=verbose)$sum/n2
  }

  if (!missing(fhat2))
  {
    fhat2 <- find.nearest.gridpts(x=rbind(x1,x2), gridx=fhat2$eval.points, f=fhat2$estimate)$fx
    psi21 <- sum(head(fhat2, n=n1))/n1
    var.fhat2 <- var(tail(fhat2, n=n2))
  }
  else
  {
    if (missing(var.fhat2))
    {
      #fhat2.x2.params <- kde.points.sum(x=x2, H=H2, eval.points=x2)
      #var.fhat2 <- (fhat2.x2.params$sumsq - fhat2.x2.params$sum^2/n2)/(n2-1)
      S2 <- var(x2)
      H2.r1 <- Hamise.mixt(mus=rep(0,d), Sigmas=S2, samp=n2, props=1, deriv.order=1)
      fhat2.r1 <- kdde(x=x2, H=H2.r1, deriv.order=1, eval.points=apply(x2, 2, mean))$estimate
      var.fhat2 <- drop(fhat2.r1 %*% S2 %*% t(fhat2.r1))
    }
    psi21 <- kde.points.sum(x=x2, H=H2, eval.points=x1, verbose=verbose)$sum/n1
  }

  ## test statistic + its parameters  
  T.hat <- drop(psi1 + psi2 - (psi12 + psi21))
  muT.hat <- (n1^(-1)*det(H1)^(-1/2) + n2^(-1)*det(H2)^(-1/2))*K0
  varT.hat <- 3*(n1*var.fhat1 + n2*var.fhat2)/(n1+n2) *(1/n1+1/n2) 
  zstat <- (T.hat-muT.hat)/sqrt(varT.hat)
  pval <- 1-pnorm(zstat)

  val <- list(Tstat=T.hat, zstat=zstat, pvalue=pval, mean=muT.hat, var=varT.hat, var.fhat1=var.fhat1, var.fhat2=var.fhat2, n1=n1, n2=n2, H1=H1, H2=H2, psi1=psi1, psi12=psi12, psi21=psi21, psi2=psi2)
  return(val)
}     
