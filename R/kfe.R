
#############################################################################
## Kernel functional estimation
#############################################################################

kfe <- function(x, G, deriv.order, inc=1, binned=FALSE, bin.par, bgridsize, deriv.vec=TRUE, add.index=TRUE, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  psir <- dmvnorm.deriv.sum(x=x, Sigma=G, deriv.order=r, inc=inc, binned=binned, bgridsize=bgridsize, deriv.vec=deriv.vec, verbose=verbose, kfe=TRUE, add.index=FALSE)
 
 if (add.index)
 {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=r, only.index=TRUE)
  
    if (deriv.vec) return(list(psir=psir, deriv.ind=ind.mat))
    else return(list(psir=psir, deriv.ind=unique(ind.mat)))
  }
  else return(psir=psir)
}


kfe.1d <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par)
{
  r <- deriv.order
  n <- length(x)
  psir <- dnorm.deriv.sum(x=x, sigma=g, deriv.order=r, inc=1, binned=binned, bin.par=bin.par, kfe=TRUE)
  if (inc==0)  psir <- (n^2*psir - n*dnorm.deriv(0, mu=0, sigma=g, deriv.order=r))/(n*(n-1))
  
  return(psir) 
}


kfe.scalar <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  if (missing(bin.par) & binned) bin.par <- binning(x=x, H=g^2*diag(d))
  
  psir <- dmvnorm.deriv.scalar.sum(x=x, sigma=g, deriv.order=r, inc=inc, kfe=TRUE, binned=binned, bin.par=bin.par, verbose=verbose)
  return(psir)
}


###############################################################################
## Plug-in unconstrained bandwidth for KFE
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


Hpi.kfe <- function(x, nstage=2, pilot, pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("currently only deriv.order=0 is implemented")
   
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  
 if (missing(pilot)) pilot <- "dscalar"
  if(!is.matrix(x)) x <- as.matrix(x)
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))  
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
  
  if (pre1=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (pre1=="sphere")
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
    if (pilot1=="dscalar")
    {
      psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
      ##h2 <- gdscalar(x=x.star, d=d, r=r-2, n=n, verbose=verbose, nstage=nstage, scv=FALSE)     
      ## gdscalar not used because it's an implicit computation without
      ## symmetriser matrices and is slower than direct computation with 
      ## symmetriser matrices. 

      A1 <- drop(t(D2K0) %*% D2K0)
      A2 <- drop(t(D2K0) %*% t(vec(diag(d)) %x% diag(d^2)) %*% psi4.ns) 
      A3 <- drop(t(psi4.ns) %*% (vec(diag(d)) %*% t(vec(diag(d))) %x% diag(d^2)) %*% psi4.ns)

      ## Special case from Chacon & Duong (2009): bias minimisation
      ##h2 <- ((4*d+8)*A1/(-d*A2 + sqrt(d^2*A2^2 + (8*d+16)*A1*A3)))^(1/(d+4))*n^(-1/(d+4))
      ## Special from Chacon & Duong (2009): bias annihilation
      h2 <- (-A1/(2*A2*n))^(1/(d+4))
      H2 <- h2^2*diag(d)   
      psi2.hat <- kfe(x=x.star, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
    else if (pilot1=="dunconstr")
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
     
 
      if (optim.fun1=="nlm")
      {
         result <- nlm(p=vech(Hstart2), f=amse2.temp, print.level=2*as.numeric(verbose))    
         H2 <- invvech(result$estimate) %*% invvech(result$estimate)
      }
      else if (optim.fun1=="optim")
      {
         result <- optim(vech(Hstart2), amse2.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
         H2 <- invvech(result$par) %*% invvech(result$par)
      }
 
      psi2.hat <- kfe(x=x, G=H2, deriv.order=r+2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
    }
  }
  else
  {
    if (pilot1=="dscalar") psi2.hat <- psins(r=r+2, Sigma=var(x.star), deriv.vec=TRUE)
    else if (pilot1=="dunconstr") psi2.hat <- psins(r=r+2, Sigma=var(x), deriv.vec=TRUE)    
  }

  if (pilot1=="dscalar") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x.star))}
  else if (pilot1=="dunconstr") {if (missing(Hstart)) Hstart <- Gns(r=r, n=n, Sigma=var(x))}
  
  ## stage 2
  amse.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    amse.val <- 1/(det(H)^(1/2)*n)*K0 + 1/2* t(vec(H)) %*% psi2.hat
    return(sum((amse.val^2))) 
  }
  
  Hstart <- matrix.sqrt(Hstart)
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=amse.temp, print.level=2*as.numeric(verbose)) 
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(vech(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }
  
  ## back-transform
  if (pilot1=="dscalar") H <- S12 %*% H %*% S12 

  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}

###############################################################################
## Plug-in diagonal bandwidth for KFE
##
## Returns
## Plug-in bandwidth
###############################################################################


 
Hpi.diag.kfe <- function(x, nstage=2, pilot, pre="scale", Hstart, binned=FALSE, bgridsize, amise=FALSE,  deriv.order=0, verbose=FALSE, optim.fun="nlm")
{
  if (deriv.order!=0) stop("currently only dervi.order=0 is implemented")

  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  K0 <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0)

  if (missing(pilot)) pilot <- "dscalar"
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))  
  pre1 <- match.arg(pre, c("scale", "sphere"))
  optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))

  if (pre1=="sphere") stop("using pre-sphering won't give diagonal bandwidth matrix")
  if (pilot1=="dunconstr") stop("unconstrained pilot selectors are not suitable for Hpi.diag.kfe")
  
  if (pre1=="scale")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }

  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=r+4, Sigma=var(x.star), deriv.vec=TRUE)
   
    if (pilot1=="dscalar")
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
  
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=diag(Hstart), f=amse.temp, print.level=2*as.numeric(verbose))    
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else if (optim.fun1=="optim")
  {
    result <- optim(diag(Hstart), amse.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }
  ## back-transform
  if (pilot1=="dscalar") H <- S12 %*% H %*% S12 
  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}











#####################################################################
## Matt Wand's version of binned kernel density derivative estimation
## Used in the feature library
##
## Computes the mth derivative of a binned
## d-variate kernel density estimate based
## on grid counts.
#############################################################

drvkde <- function(x,drv,bandwidth,gridsize,range.x,binned=FALSE,se=TRUE, w)
{  
   d <- length(drv)
   if (d==1) x <- as.matrix(x)

   ## Rename common variables
   h <- bandwidth
   tau <- 4 + max(drv)    
   if (length(h)==1) h <- rep(h,d)

   if (missing(gridsize))
     if (!binned)   ## changes 16/02/2009
     {  
       if (d==1) gridsize <- 401
       else if (d==2) gridsize <- rep(151,d)
       else if (d==3) gridsize <- rep(51, d)
       else if (d==4) gridsize <- rep(21, d)
     }
     else
     {
       if (d==1) gridsize <- dim(x)[1]
       else gridsize <- dim(x)
     }

   if(missing(w)) w <- rep(1,nrow(x))
   ## Bin the data if not already binned
  
   if (missing(range.x)) 
   {
     range.x <- list()
     for (id in 1:d)
       range.x[[id]] <- c(min(x[,id])-tau*h[id],max(x[,id])+tau*h[id])  
   }
   
   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))
   
   M <- gridsize
   gpoints <- list()

   for (id in 1:d)
     gpoints[[id]] <- seq(a[id],b[id],length=M[id])

   if (binned==FALSE)
   {
     if (d==1) gcounts <- linbin.ks(x,gpoints[[1]], w=w)
     if (d==2) gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]], w=w)
     if (d==3) gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]], w=w)
     if (d==4) gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]], w=w)
   }
   else
     gcounts <- x

   n <- sum(gcounts)

   kapmid <- list()
   for (id in (1:d))
   {
     ## changes to Lid 13/02/2009
     Lid <- max(min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),M[id]),d)
     lvecid <- (0:Lid)
     facid  <- (b[id]-a[id])/(h[id]*(M[id]-1))
     argid <- lvecid*facid
     kapmid[[id]] <- dnorm(argid)/(h[id]^(drv[id]+1))
     hmold0 <- 1
     hmold1 <- argid
     if (drv[id]==0) hmnew <- 1
     if (drv[id]==1) hmnew <- argid
     if (drv[id] >= 2) 
       for (ihm in (2:drv[id])) 
       {
         hmnew <- argid*hmold1 - (ihm-1)*hmold0
         hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
         hmold1 <- hmnew    # by recurrence.
       }
     kapmid[[id]] <- hmnew*kapmid[[id]]*(-1)^drv[id]
   }
  
   if (d==1) kappam <- kapmid[[1]]/n
   if (d==2) kappam <- outer(kapmid[[1]],kapmid[[2]])/n
   if (d==3) kappam <- outer(kapmid[[1]],outer(kapmid[[2]],kapmid[[3]]))/n
   if (d==4) kappam <- outer(kapmid[[1]],outer(kapmid[[2]],outer(kapmid[[3]],kapmid[[4]])))/n
  
   if (!any(c(d==1,d==2,d==3,d==4))) stop("only for d=1,2,3,4")

   if (d==1) 
   {
     est <- symconv.ks(kappam,gcounts,skewflag=(-1)^drv)
     if (se) est.var <- ((symconv.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }

   if (d==2) 
   {     
     est <- symconv2D.ks(kappam,gcounts,skewflag=(-1)^drv)
     if (se) est.var <- ((symconv2D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
   }
     
   if (d==3)
   {
     est <- symconv3D.ks(kappam,gcounts,skewflag=(-1)^drv) 
     if (se) est.var <- ((symconv3D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)
   }
     
   if (d==4)
   {
     est <- symconv4D.ks(kappam,gcounts,skewflag=(-1)^drv) 
     if (se) est.var <- ((symconv4D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }
   
   if (se)
   {
     est.var[est.var<0] <- 0
     return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
   }
   else if (!se)
     return(list(x.grid=gpoints,est=est))
}

