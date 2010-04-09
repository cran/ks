###############################################################################
### Multivariate kernel density derivative estimate 
###############################################################################

kdde <- function(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w)
{
  r <- deriv.order

  if (is.vector(x))
  {
    if (missing(H)) {d <- 1; n <- length(x)}
    else
    {
      if (is.vector(H)) { d <- 1; n <- length(x)}
      else {x <- matrix(x, nrow=1); d <- ncol(x); n <- nrow(x)}
    }
  }
  else {d <- ncol(x); n <- nrow(x)}

  if (!missing(w))
    if (!(identical(all.equal(sum(w), n), TRUE)))
    {
      warning("Weights don't sum to sample size - they have been scaled accordingly\n")
      w <- w*n/sum(w)
    }

  if (missing(w)) w <- rep(1,n)


  ## compute binned estimator
  if (binned)
  {
    if (!missing(eval.points))
      stop("Both binned=TRUE and eval.points are non-empty.")
    
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
  
    if (positive & is.vector(x))
    {
      y <- log(x)
      fhat <- kdde.binned(x=y, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, single.deriv=FALSE)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
      fhat <- kdde.binned(x=x, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, single.deriv=FALSE)
  }
  else
  {
    ## compute exact (non-binned) estimator
    if (missing(gridsize)) gridsize <- default.gridsize(d)
    
    ## 1-dimensional    
    if (d==1)
    {
      if (!missing(H) & !missing(h))
        stop("Both H and h are both specified.")
      
      if (missing(h))
        h <- sqrt(H)

      if (missing(eval.points))
        fhat <- kdde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r)
      else
        fhat <- kdde.points.1d(x=x, h=h, eval.points=eval.points, w=w, deriv.order=r)
    }
    ## multi-dimensional
    else
    {   
      if (is.data.frame(x)) x <- as.matrix(x)
      
      if (missing(eval.points))
      {
        if (d==2) ##stop("Not yet implemented for 2 dimensions")
          fhat <- kdde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r)
        else if (d == 3) stop("Not yet implemented for 3 dimensions")
          #fhat <- kde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w) 
        else 
          stop("Need to specify eval.points for more than 3 dimensions")
      }
      else
        fhat <- kdde.points(x, H, eval.points, w=w, deriv.order=r)
    }

  }
  fhat$binned <- binned

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

  ## rearrange list fields
 
      
  class(fhat) <- "kdde"

  return(fhat)
 }



###############################################################################
## Multivariate binned kernel density derivative estimate
## for single partial derivative ONLY
###############################################################################

### Diagonal H only
kdde.binned <- function(x, H, h, deriv.order, bgridsize, xmin, xmax, bin.par, w, single.deriv=TRUE)
{
  r <- deriv.order
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
    bin.par <- binning(x=x, H=H, h, bgridsize, xmin, xmax, supp=3.7+max(r), w=w)
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

  if (d==1)
  {  
    fhat.grid <- drvkde(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(H)), binned=TRUE, range.x=range.x, se=FALSE, gridsize=bgridsize)
    eval.points <- fhat.grid$x.grid
    fhat.grid <- fhat.grid$est
  }
  else
  {
    if (single.deriv)
    {
      fhat.grid <- drvkde(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(H)), binned=TRUE, range.x=range.x, se=FALSE, gridsize=bgridsize)
      eval.points <- fhat.grid$x.grid
      fhat.grid <- fhat.grid$est
      ind.mat <- r
    }
    else
    {  
      ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r,  index.only=TRUE)
      fhat.grid <- list()
      
      ## Needs to be optimised here to avoid duplicated computation of density derviatives
      for (r2 in 1:nrow(ind.mat))
      {
        fhat.gridr2 <- drvkde(x=bin.par$counts, drv=ind.mat[r2,], bandwidth=sqrt(diag(H)), binned=TRUE, range.x=range.x, se=FALSE, gridsize=bgridsize)
        fhat.grid[[r2]] <- fhat.gridr2$est
      }
      eval.points <- fhat.grid$x.grid
    }
  }
  
  if (missing(x)) x <- NULL
  
  if (d==1)
    fhat <- list(x=x, eval.points=unlist(eval.points), estimate=fhat.grid, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=r)
  else
    fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat)

  class(fhat) <- "kdde"
  
  return(fhat)
}




##############################################################################################
#### Univariate kernel density derivative estimate on a grid
##############################################################################################

kdde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype, w, deriv.order=0)
{
  r <- deriv.order
  if (r==0)
    fhatr <- kde(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, adj.positive=adj.positive, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w)
  else
  {  
    if (missing(xmin)) xmin <- min(x) - h*supp
    if (missing(xmax)) xmax <- max(x) + h*supp
    if (missing(gridtype)) gridtype <- "linear"
  
    y <- x
    gridtype1 <- tolower(substr(gridtype,1,1))
    if (gridtype1=="l")
    {
      gridy <- seq(xmin, xmax, length=gridsize)
      gridtype.vec <- "linear"
    }
    else if (gridtype1=="s")
    {
      gridy.temp <- seq(sign(xmin)*sqrt(abs(xmin)), sign(xmax)*sqrt(abs(xmax)), length=gridsize)
      gridy <- sign(gridy.temp) * gridy.temp^2
      gridtype.vec <- "sqrt"
    }
  
    n <- length(y)
    est <- dnorm.deriv.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=w/n, deriv.order=r)
    fhatr <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=deriv.order)
      
    class(fhatr) <- "kde"
  }
  
  return(fhatr)
}



##############################################################################################
## Bivariate kernel density derivative estimate on a grid
## Computes all mixed partial derivatives for a given deriv.order
##############################################################################################

kdde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0)
{
  d <- 2
  r <- deriv.order
  if (r==0)
    fhatr <- kde(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w)
  else
  {  
    ## initialise grid 
    n <- nrow(x)
    if (is.null(gridx))
      gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
    
    suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
    
    if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    

    nderiv <- d^r
    fhat.grid <- list()
    for (k in 1:nderiv)
      fhat.grid[[k]] <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))

    S2r <- Sdr(d=d, r=r)
    for (i in 1:n)
    {
      ## compute evaluation points 
      eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
      eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
      eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
      eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
      eval.x.len <- length(eval.x)
      eval.pts <- permute(list(eval.x, eval.y))

      ## Create list of matrices for different partial derivatives
      fhat <- dmvnorm.deriv(x=eval.pts, mu=x[i,], Sigma=H, deriv.order=r, Sdr=S2r)
      
      ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
      for (k in 1:nderiv)
        for (j in 1:length(eval.y))
          fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] <- fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] + w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len),k]
    }
    
    for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
    gridx1 <- list(gridx[[1]], gridx[[2]]) 

    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, Sdr=S2r, index.only=TRUE)
    fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
  }

  return(fhatr)
}


#################################################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
#################################################################################################

kdde.points <- function(x, H, eval.points, w, deriv.order=0) 
{
  n <- nrow(x)
  Hs <- numeric(0)
  for (i in 1:n)
    Hs <- rbind(Hs, H)
  r <- deriv.order
  d <- ncol(H)
  fhat <- dmvnorm.deriv.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n, deriv.order=r)
  
  ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, Sdr=Sdr(d=d, r=r), index.only=TRUE)
  return(list(x=x, eval.points=eval.points, estimate=fhat, H=H, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat))
}

kdde.points.1d <- function(x, h, eval.points, w, deriv.order=0) 
{
  r <- deriv.order
  n <- length(x)
  fhat <- dnorm.deriv.mixt(x=eval.points, mus=x, sigmas=rep(h,n), props=w/n, deriv.order=r)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=r))
}

#####################################################################
### Matt Wand's version of binned kernel density derivative estimation
###
### Computes the mth derivative of a binned
### d-variate kernel density estimate based
### on grid counts.
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
    psir <- dnorm.deriv.sum(x=x, sigma=g, deriv.order=r, inc=inc, binned=FALSE, kfe=TRUE)
  else
    #psir <- dnorm.deriv.sum(bin.par=bin.par, sigma=g, r=r, inc=inc, kfe=TRUE, binned=TRUE)
  { psir <- bkfe(x=bin.par$counts, bandwidth=g, drv=r, binned=TRUE, range.x=range(bin.par$eval.points)); if (inc==0)  psir <- psir - dnorm.deriv(0,mu=0,sigma=g, deriv.order=r)/sum(bin.par$counts)} 
  
  return(psir) 
}

kfe <- function(x, G, r, inc=1, binned=FALSE, bin.par, diff=FALSE)
{
  if (!binned)
    psir <- dmvnorm.deriv.sum(x=x, Sigma=G, deriv.order=r, inc=inc, diff=diff, kfe=TRUE, binned=FALSE)
  else
    psir <- dmvnorm.deriv.sum(bin.par=bin.par, Sigma=G, deriv.order=r, inc=inc, kfe=TRUE, binned=TRUE) 

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
      return(mat.Kpow(args,6)-15*mat.Kprod(mat.Kpow(args,4),rep(1,nrow(args))%x%t(vId)) + 45*mat.Kprod(mat.Kpow(args,2),rep(1,nrow(args))%x%t(Kpow(vId,2)))-15*rep(1,nrow(args))%x%t(Kpow(vId,3)))
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

    vecPsi6 <- (-1)^6*det(Grinv12)*Kpow(Grinv12,6)%*%vecPsi6

    return (vecPsi6)
  }

  
  if (r==4)
  {
    vecPsi4 <- rep(0,d^4)
    K4 <- function(args)
    {
      return(mat.Kpow(args,4)-6*mat.Kprod(mat.Kpow(args,2),rep(1,nrow(args))%x%t(vId))+3*rep(1,nrow(args))%x%t(Kpow(vId,2)))
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

    vecPsi4<-(-1)^4*det(Grinv12)*Kpow(Grinv12,4)%*%vecPsi4

    return(vecPsi4)
  }
}



