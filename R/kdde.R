###############################################################################
### Multivariate kernel density derivative estimate 
###############################################################################

kdde <- function(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, deriv.vec=TRUE, verbose=FALSE)
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
    if (!missing(eval.points)) stop("Both binned=TRUE and eval.points are non-empty.")
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    ##if (!identical(diag(diag(H)), H) & d > 1) stop("Binned estimation defined for diagonal H only")
    
    if (positive & is.vector(x))
    {
      y <- log(x)
      fhat <- kdde.binned(x=y, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
      fhat <- kdde.binned(x=x, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.vec=deriv.vec, verbose=verbose)
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
        if (d==2) 
          fhat <- kdde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r, deriv.vec=deriv.vec)
        else if (d == 3) stop("Not yet implemented for 3 dimensions")
          ##fhat <- kde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w) 
        else 
          stop("Need to specify eval.points for more than 3 dimensions")
      }
      else
        fhat <- kdde.points(x=x, H=H, eval.points=eval.points, w=w, deriv.order=r, deriv.vec=deriv.vec)
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
  class(fhat) <- "kdde"

  return(fhat)
 }



###############################################################################
## Multivariate binned kernel density derivative estimate
###############################################################################

kdde.binned <- function(x, H, h, deriv.order, bgridsize, xmin, xmax, bin.par, w, deriv.vec=TRUE, deriv.index, Sdr.mat, verbose=FALSE)
{
  r <- deriv.order
  if (length(r)>1) stop("deriv.order should be a non-negative integer.")
  
  ## linear binning
  if (missing(bin.par))
  {
    if (is.vector(x)) {d <- 1; n <- length(w)}
    else {d <- ncol(x); n <- nrow(x)}

    if (missing(w)) w <- rep(1,n)
    
    if (d==1)
      if (missing(H)) { H <- as.matrix(h^2)} 
      else {h <- sqrt(H); H <- as.matrix(H)}

    if (d==1) Hd <- H else Hd <- diag(diag(H))
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    
    bin.par <- binning(x=x, H=Hd, h=h, bgridsize, xmin, xmax, supp=3.7+max(r), w=w)
  }
  else
  {
    if (!is.list(bin.par$eval.points)) { d <- 1; bgridsize <- length(bin.par$eval.points)}
    else  { d <- length(bin.par$eval.points); bgridsize <- sapply(bin.par$eval.points, length)} 

    w <- bin.par$w
    if (d==1)
      if (missing(H)) H <- as.matrix(h^2)
      else {h <- sqrt(H); H <- as.matrix(H)}
  }

 
  if (d==1)
  {
    fhat <- kdde.binned.1d(h=h, deriv.order=r, bin.par=bin.par)
    eval.points <- fhat$eval.points
    est <- fhat$estimate
  }
  else
  {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE, deriv.vec=TRUE)
    fhat.grid <- kdde.binned.nd(H=H, deriv.order=r, bin.par=bin.par, Sdr.mat=Sdr.mat, verbose=verbose)
  }

  if (missing(x)) x <- NULL
  
  if (d==1)
    fhat <- list(x=x, eval.points=unlist(eval.points), estimate=est, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=r)
  else
  {
    if (r==0)
      fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate[[1]], H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat)
    else
      fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate, H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat)
  }
  
  class(fhat) <- "kdde"
  
  return(fhat)
}

kdde.binned.1d <- function(h, deriv.order, bin.par)
{
  r <- deriv.order
  n <- sum(bin.par$counts)
  a <- min(bin.par$eval.points)
  b <- max(bin.par$eval.points)
  M <- length(bin.par$eval.points)
  L <- min(ceiling((4+r)*h*(M-1)/(b-a)), M-1)
  Keval <- dnorm.deriv(x=(b-a)*(0:L)/(M-1), mu=0, sigma=h, deriv.order=r)/n
  est <- symconv.ks(Keval, bin.par$counts, skewflag=(-1)^r)

  return(list(eval.points=bin.par$eval.points, estimate=est))
}

kdde.binned.nd <- function(H, deriv.order, bin.par, Sdr.mat, verbose=FALSE)
{
  d <- ncol(H)
  r <- deriv.order
  n <- sum(bin.par$counts)
  a <- sapply(bin.par$eval.points, min)
  b <- sapply(bin.par$eval.points, max)
  M <- sapply(bin.par$eval.points, length)
  L <- pmin(ceiling((4+r)*max(sqrt(abs(diag(H))))*(M-1)/(b-a)), M-1)

  if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d, r=r)
  if (d==2) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2])
  if (d==3) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2], (b[3]-a[3])*(0:L[3])/M[3])
  if (d==4) xgrid <- expand.grid((b[1]-a[1])*(0:L[1])/M[1], (b[2]-a[2])*(0:L[2])/M[2], (b[3]-a[3])*(0:L[3])/M[3], (b[4]-a[4])*(0:L[4])/M[4])
  ##nxgrid <- nrow(xgrid)

  deriv.index <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, Sdr.mat=Sdr.mat, only.index=TRUE) ##Keval$deriv.ind
  deriv.index.minimal <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, Sdr.mat=Sdr.mat, only.index=TRUE, deriv.vec=FALSE)

  Keval <- dmvnorm.deriv(x=xgrid, mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, Sdr.mat=Sdr.mat, deriv.vec=FALSE)
  Keval <- Keval$deriv/n
  if (r==0) Keval <- as.matrix(Keval, ncol=1)
  
  est <- list()
  if (verbose) pb <- txtProgressBar() 

  ## loop over only unique partial derivative indices
  nderiv <- nrow(deriv.index.minimal)
  if (!(is.null(nderiv)))
    for (s in 1:nderiv)
    {
      deriv.rep.index <- which.mat(deriv.index.minimal[s,], deriv.index)
      Kevals <- array(Keval[,s], dim=L+1)
      if (r==0) sf <- rep(1,d)
      else sf <- (-1)^deriv.index.minimal[s,]
      if (d==2) est.temp <- symconv2D.ks(Kevals, bin.par$counts, skewflag=sf)
      if (d==3) est.temp <- symconv3D.ks(Kevals, bin.par$counts, skewflag=sf)
      if (d==4) est.temp <- symconv4D.ks(Kevals, bin.par$counts, skewflag=sf)
      for (s2 in 1:length(deriv.rep.index)) est[[deriv.rep.index[s2]]] <-  zapsmall(est.temp)
      if (verbose) setTxtProgressBar(pb, s/nderiv)
    }
  else
  {
    for (s in 1:ncol(Keval))
    {  
      Kevals <- array(Keval[,s], dim=L+1)
      if (r==0) sf <- rep(1,d)
      else sf <- (-1)^deriv.index[s,]
      if (d==2) est[[s]] <- symconv2D.ks(Kevals, bin.par$counts, skewflag=sf)
      if (d==3) est[[s]] <- symconv3D.ks(Kevals, bin.par$counts, skewflag=sf)
      if (d==4) est[[s]] <- symconv4D.ks(Kevals, bin.par$counts, skewflag=sf)
      est[[s]] <- zapsmall(est[[s]])
      if (verbose) setTxtProgressBar(pb, s/ncol(Keval))
    }
  }
  if (verbose) close(pb)
  

  return(list(eval.points=bin.par$eval.points, estimate=est, deriv.order=r))
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

kdde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0, deriv.vec=TRUE)
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
      fhat <- dmvnorm.deriv(x=eval.pts, mu=x[i,], Sigma=H, deriv.order=r, Sdr.mat=S2r)
      
      ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
      for (k in 1:nderiv)
        for (j in 1:length(eval.y))
          fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] <- fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] + w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len),k]
    }
    
    for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
    gridx1 <- list(gridx[[1]], gridx[[2]]) 

    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, Sdr.mat=S2r, only.index=TRUE)

    if (!deriv.vec)
    {
      fhat.grid.vech <- list()
      deriv.ind <- unique(ind.mat)

      for (i in 1:nrow(deriv.ind))
      {
        which.deriv <- which.mat(deriv.ind[i,], ind.mat)[1]
        fhat.grid.vech[[i]] <- fhat.grid[[which.deriv]]
      }
      ind.mat <- deriv.ind
      fhat.grid <- fhat.grid.vech
    }

    fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
  }

  return(fhatr)
}


#################################################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
#################################################################################################

kdde.points.1d <- function(x, h, eval.points, w, deriv.order=0) 
{
  r <- deriv.order
  n <- length(x)
  fhat <- dnorm.deriv.mixt(x=eval.points, mus=x, sigmas=rep(h,n), props=w/n, deriv.order=r)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=r))
}


kdde.points <- function(x, H, eval.points, w, deriv.order=0, deriv.vec=TRUE) 
{
  n <- nrow(x)
  Hs <- numeric(0)
  for (i in 1:n)
    Hs <- rbind(Hs, H)
  r <- deriv.order
  fhat <- dmvnorm.deriv.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n, deriv.order=r, deriv.vec=deriv.vec, add.index=TRUE)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat$deriv, H=H, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=fhat$deriv.ind))
}



#############################################################################
## Kernel functional estimation
#############################################################################

kfe.1d <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par)
{
  r <- deriv.order
  psir <- dnorm.deriv.sum(x=x, sigma=g, deriv.order=r, inc=inc, binned=binned, bin.par=bin.par, kfe=TRUE)
  if (inc==0)  psir <- psir - dnorm.deriv(0, mu=0, sigma=g, deriv.order=r)/length(x)
  
  return(psir) 
}

kfe <- function(x, G, deriv.order, inc=1, binned=FALSE, bin.par, bgridsize, double.loop=FALSE, deriv.vec=TRUE, add.index=TRUE, Sdr.mat, verbose=FALSE, kfold=1)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)

  ##if (missing(double.loop)) double.loop <- (n*d^r/kfold > 1e6)
  if (kfold > 1)
  {
    kfold.group <- sample(size=n, 1:kfold, replace=TRUE)
    psir <- 0
    for (k in 1:kfold)
      psir <- psir + dmvnorm.deriv.sum(x=x[which(kfold.group==k),], Sigma=G, deriv.order=r, inc=inc, binned=FALSE, deriv.vec=deriv.vec, add.index=FALSE, double.loop=double.loop, Sdr.mat=Sdr.mat, verbose=verbose, kfe=TRUE)
    psir <- psir/kfold
  }
  else
  {  
    psir <- dmvnorm.deriv.sum(x=x, Sigma=G, deriv.order=r, inc=inc, binned=binned, double.loop=double.loop, bin.par=bin.par, bgridsize=bgridsize, deriv.vec=deriv.vec, verbose=verbose, Sdr.mat=Sdr.mat, kfe=TRUE)
  }
  psir <- drop(psir)
  
  
  if (add.index)
  {
    ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=r, only.index=TRUE)
  
    if (deriv.vec) return(list(psir=psir, deriv.ind=ind.mat))
    else return(list(psir=psir, deriv.ind=unique(ind.mat)))
  }
  else return(psir=psir)
}


kfe.scalar <- function(x, g, deriv.order, inc=1, binned=FALSE, bin.par, double.loop=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  if (missing(bin.par) & binned) bin.par <- binning(x=x, H=g^2*diag(d))
  
  psir <- dmvnorm.deriv.scalar.sum(x=x, sigma=g, deriv.order=r, inc=inc, kfe=TRUE, binned=binned, double.loop=double.loop, bin.par=bin.par)
  return(psir)
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

