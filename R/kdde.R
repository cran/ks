###############################################################################
### Multivariate kernel density derivative estimate 
###############################################################################

kdde <- function(x, H, h, r, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=TRUE, bgridsize, positive=FALSE, adj.positive)
{
  if (is.vector(x))
  {
    if (missing(H)) d <- 1
    else
    {
      if (is.vector(H)) d <- 1
      else {x <- matrix(x, nrow=1); d <- ncol(x)}
    }
  }
  else d <- ncol(x)
  
  ## compute binned estimator
  if (binned)
  {
    if (!missing(eval.points))
      stop("Both binned=TRUE and eval.points are non-empty.")
    
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
  
    if (positive & is.vector(x))
    {
      y <- log(x)
      fhat <- kdde.binned(x=y, H=H, h=h, r=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
      fhat <- kdde.binned(x=x, H=H, h=h, r=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax)
  }
  else
  {
    ## compute exact (non-binned) estimator
    stop("Non-binned estimators for density derivatives not yet implemented")
  }
  
  fhat$binned <- binned
  ##fhat$gridtype <- gridtype

  ## add variable names
  if (is.vector(x))
  {
    d <- 1
    x.names <- deparse(substitute(x))
  }
  else
  {  
    d <- ncol(x)
    x.names <- colnames(x)
    if (is.null(x.names))
    {
      x.names <- strsplit(deparse(substitute(x)), "\\[")[[1]][1]
      x.names <- paste(x.names, "[,", 1:d,"]",sep="") 
    }
  }
  fhat$names <- x.names
  class(fhat) <- "kde"
  

  return(fhat)
 }



###############################################################################
### Multivariate binned kernel density derivative estimate 
###############################################################################

### Diagonal H only
kdde.binned <- function(x, H, h, r, bgridsize, xmin, xmax, bin.par)
{
  ## linear binning
  
  if (missing(bin.par))
  {
    if (is.vector(x)) d <- 1
    else d <- ncol(x)

    if (d==1)
      if (missing(H)) H <- as.matrix(h^2)
      else {h <- sqrt(H); H <- as.matrix(H)}

    if (!identical(diag(diag(H)), H) & d > 1)
      stop("Binned estimation defined for diagonal Sigma only")
  
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    bin.par <- binning(x=x, H=H, h, bgridsize, xmin, xmax, supp=3.7+max(r))
  }
  else
  {
    if (!is.list(bin.par$eval.points)) { d <- 1; bgridsize <- length(bin.par$eval.points)}
    else  { d <- length(bin.par$eval.points); bgridsize <- sapply(bin.par$eval.points, length)} 
    
    if (d==1)
      if (missing(H)) H <- as.matrix(h^2)
      else {h <- sqrt(H); H <- as.matrix(H)}
  }
  
  if (d==1)
    range.x <- list(range(bin.par$eval.points))
  else
    range.x <- lapply(bin.par$eval.points, range)

 
  fhat.grid <- drvkde(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(H)), binned=TRUE, range.x=range.x, se=FALSE, gridsize=bgridsize)
  eval.points <- fhat.grid$x.grid
  fhat.grid <- fhat.grid$est
  
  if (missing(x)) x <- NULL
  
  if (d==1)
    fhat <- list(x=x, eval.points=unlist(eval.points), estimate=fhat.grid, H=h^2, h=h)
  else
    fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, H=H)
  
  return(fhat)
}


#####################################################################
### Matt Wand's version of binned kernel density derivative estimation
###
### Computes the mth derivative of a binned
### d-variate kernel density estimate based
### on grid counts.
#############################################################

drvkde <- function(x,drv,bandwidth,gridsize,range.x,binned=FALSE,se=TRUE)
{  
   d <- length(drv)
   if (d==1) x <- as.matrix(x)

   ## Rename common variables
   h <- bandwidth
   tau <- 4 + max(drv)    
   if (length(h)==1) h <- rep(h,d)

   if (missing(gridsize))
     if (d==1) gridsize <- 401
     else if (d==2) gridsize <- rep(151,d)
     else if (d==3) gridsize <- rep(51, d)
     else if (d==4) gridsize <- rep(21, d)
   
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
     if (d==1) gcounts <- linbin.ks(x,gpoints[[1]])
     if (d==2) gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]])
     if (d==3) gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]])
     if (d==4) gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]])
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
     ##kappam <- as.vector(kappam)
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


#############################################################################
## Kernel functional estimation
#############################################################################

kfe.1d <- function(x, g, r, inc=1, binned=FALSE, bin.par)
{
  if (!binned)
    psir <- dnorm.deriv.sum(x=x, sigma=g, r=r, inc=inc, binned=FALSE, kfe=TRUE)
  else
    #psir <- dnorm.deriv.sum(bin.par=bin.par, sigma=g, r=r, inc=inc, kfe=TRUE, binned=TRUE)
  { psir <- bkfe(x=bin.par$counts, bandwidth=g, drv=r, binned=TRUE, range.x=range(bin.par$eval.points)); if (inc==0)  psir <- psir - dnorm.deriv(0,mu=0,sigma=g, r=r)/sum(bin.par$counts)} 
  
  return(psir) 
}

kfe <- function(x, G, r, inc=1, binned=FALSE, bin.par, diff=FALSE)
{ 
  if (!binned)
    psir <- dmvnorm.deriv.sum(x=x, Sigma=G, r=r, inc=inc, diff=diff, kfe=TRUE, binned=FALSE)
  else
    psir <- dmvnorm.deriv.sum(bin.par=bin.par, Sigma=G, r=r, inc=inc, kfe=TRUE, binned=TRUE) 

  return(psir)
}



#############################################################################
## Estimation of psi_r (no binning, arbitrary d) for G = g^2 I
##
## Code by Jose E. Chacon. Received  04/09/2007
## x - data matrix (or differences matrix)
## r - derivative
## g - scalar pilot bandwidth
## diff - flag for x is data or differences matrix
## upper - compute upper diagonal of differences matrix
#############################################################################
    
##psir.hat  function(x, r, g, diff=FALSE, upper=diff)
kfe.scalar <- function(x, g, r, diff=FALSE)
{    
  if(!diff)
  {
    n <- nrow(x)
    d <- ncol(x)
    if (d != length(r))
      stop("The length of r must equal the number of columns of x")
    
    difs <- differences(x, upper=FALSE)}
  else
  {
    n <- sqrt(nrow(x)) ##(-1 + sqrt(1+8*nrow(x)))/2
    d <- ncol(x)
    difs <- x
  }
  
  sderiv <- sum(r)
  arg <- difs/g
  darg <- dmvnorm(arg)/(g^(sderiv+d))
  for (j in 1:d)
  {
    hmold0 <- 1
    hmold1 <- arg[,j]
    hmnew <- 1
    if (r[j] ==1){hmnew<-hmold1}
    if (r[j] >= 2) ## Multiply by the corresponding Hermite polynomial, coordinate-wise, using Fact C.1.4 in W&J (1995) and Willink (2005, p.273)
      for (i in (2:r[j]))
      {
        hmnew <- arg[,j] * hmold1 - (i - 1) * hmold0
        hmold0 <- hmold1
        hmold1 <- hmnew
      }
    darg <- hmnew * darg
  }

  psir <- mean(darg)*(-1)^sderiv      
  return(psir)
}


#############################################################################
## Estimation of vec Psi_r (no binning, arbitrary d) for unconstrained G
#############################################################################

vecPsir <- function(x, Gr, Sdr, r, upper, nlim=1e4)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  ## matrix of differences - upper triangular form is available
  difs <- differences(x, upper=upper)
  
  Id1 <- diag(d)
  vId <- vec(Id1)
  Gr12 <- matrix.sqrt(Gr)
  Grinv12 <- chol2inv(chol(Gr12))

  args <- t(Grinv12%*%t(difs))

  if (r==6)
  {  
    vecPsi6 <- rep(0,d^6)

    K6 <- function(args)
    {
      return(mat.K.pow(args,6)-15*mat.Kprod(mat.K.pow(args,4),rep(1,nrow(args))%x%t(vId)) + 45*mat.Kprod(mat.K.pow(args,2),rep(1,nrow(args))%x%t(K.pow(vId,2)))-15*rep(1,nrow(args))%x%t(K.pow(vId,3)))
    }

    ## split into blocks because of memory limitations
    if (nrow(args) > nlim)
    {
      num.loops <- nrow(args) %/% nlim
      for (j in 1:num.loops)
      {
        args.temp <- args[((j-1)* nlim+1):(j*nlim),]
        hmnew.temp <- K6(args.temp)
        vecPsi6  <- vecPsi6 + colSums(dmvnorm(args.temp)*hmnew.temp)
        ##cat(j, " ")
      }
    
      if (nrow(args) %% nlim >0)
      {
        args.temp <- args[(num.loops*nlim+1):nrow(args),]
        hmnew.temp <- K6(args.temp)
        vecPsi6  <- vecPsi6 + colSums(dmvnorm(args.temp)*hmnew.temp)
      }
    }
    else
    {
      hmnew <- K6(args)
      vecPsi6 <- colSums(dmvnorm(args)*hmnew)
    }
    ##cat("\n")

    if (upper)
    {
      ## adjust for upper triangular form of differences matrix
      args0 <- t(rep(0,d))
      hmnew0 <- K6(args0)
      vecPsi6 <- 2*vecPsi6 - as.vector(n*dmvnorm(args0)*hmnew0)
      vecPsi6 <- Sdr%*%vecPsi6/n^2
    }
    else
      vecPsi6 <- Sdr%*%vecPsi6/nrow(args)

    vecPsi6 <- (-1)^6*det(Grinv12)*K.pow(Grinv12,6)%*%vecPsi6

    return (vecPsi6)
  }

  
  if (r==4)
  {
    vecPsi4 <- rep(0,d^4)
    K4 <- function(args)
    {
      return(mat.K.pow(args,4)-6*mat.Kprod(mat.K.pow(args,2),rep(1,nrow(args))%x%t(vId))+3*rep(1,nrow(args))%x%t(K.pow(vId,2)))
    }
    
    if (nrow(args)>nlim)
    {
      num.loops <- nrow(args) %/% nlim
      for (j in 1:num.loops)
      {
        args.temp <- args[((j-1)* nlim+1):(j*nlim),]
        hmnew.temp <- K4(args.temp)
        vecPsi4  <- vecPsi4 + colSums(dmvnorm(args.temp)*hmnew.temp)
      }
    
      if (nrow(args) %% nlim >0)
      {
        args.temp <- args[(num.loops*nlim+1):nrow(args),]
        hmnew.temp <- K4(args.temp)
        vecPsi4  <- vecPsi4 + colSums(dmvnorm(args.temp)*hmnew.temp)
      }   
    }
    else
    {
      hmnew <- K4(args)
      vecPsi4 <- colSums(dmvnorm(args)*hmnew)
    }
    
    if (upper)
    {
      args0 <- t(rep(0,d))
      hmnew0 <- K4(args0)
      vecPsi4 <- 2*vecPsi4 - as.vector(n*dmvnorm(args0)*hmnew0)
      vecPsi4 <- Sdr%*%vecPsi4/n^2
    }
    else
      vecPsi4 <- Sdr%*%vecPsi4/nrow(args)

    vecPsi4<-(-1)^4*det(Grinv12)*K.pow(Grinv12,4)%*%vecPsi4

    return(vecPsi4)
  }
}



