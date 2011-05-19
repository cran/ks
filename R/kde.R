##################################################################################
## Multivariate kernel density estimators
#################################################################################


##################################################################################
## Generate grid over a set of points
##
## Parameters
## x - data points
## H - bandwidth matrix
## tol - tolerance = extra coverage exceeding the range of x   
## gridsize - number of points for each direction
##
## Returns
## gridx - list of intervals, one for each co-ord direction so that
##         gridx[[1]] x gridx[[2]] x ... x gridx[[d]] is the grid
## stepsize - vector of step sizes 
##################################################################################

make.grid.ks <- function(x, H, tol, gridsize, xmin, xmax, gridtype)
{
  d <- ncol(x)
  tol.H <-  tol * diag(H)
  if (missing(xmin))
     xmin <- apply(x, 2, min) - tol.H
  if (missing(xmax))
     xmax <- apply(x, 2, max) + tol.H
  
  stepsize <- rep(0, d)
  gridx <- numeric(0)
  if (length(gridsize)==1)
    gridsize <- rep(gridsize, d)

  
  if (missing(gridtype))
   gridtype <- rep("linear", d)

  gridtype.vec<- rep("", d)
  
  for (i in 1:d)
  {
    gridtype1 <- tolower(substr(gridtype[i],1,1))
    if (gridtype1=="l")
    {  
      gridx <- c(gridx, list(seq(xmin[i], xmax[i], length=gridsize[i])))
      stepsize[i] <- abs(gridx[[i]][1] - gridx[[i]][2])
      gridtype.vec[i] <- "linear"
    }
    else if (gridtype1=="s")
    {
      gridx.temp <- seq(sign(xmin[i])*sqrt(abs(xmin[i])), sign(xmax[i])*sqrt(abs(xmax[i])), length=gridsize[i])
      gridx <- c(gridx, list(sign(gridx.temp) * gridx.temp^2))
      stepsize[i] <- NA
      gridtype.vec[i] <- "sqrt"
    }
  }
  
  gridx <- c(gridx, list(stepsize = stepsize, gridtype=gridtype.vec))
    
  return(gridx)
}  


######################################################################################
## Generate kernel (rectangular) support at data point
## 
## Parameters
## x - data points
## H - bandwidth matrix
## tol - tolerance = extra coverage exceeding the range of x 
##
## Returns
## list of min and max points of support (here we parameterise rectangles
## by their min = lower left co-ord and max = upper right coord)
#####################################################################################

make.supp <- function(x, H, tol)
{
  n <- nrow(x)
  d <- ncol(x)
  tol.H <- tol * diag(H)
  xmin <- matrix(0, nrow=n, ncol=d)
  xmax <- matrix(0, nrow=n, ncol=d)

  for (i in 1:n)
  {
    xmin[i,] <- x[i,] - tol.H
    xmax[i,] <- x[i,] + tol.H 
  }
           
  return(list(xmin = xmin, xmax = xmax))
}


######################################################################################
## Find the grid points contained in kernel support rectangles 
##
## Parameters
## gridx - grid (list of subdivided intervals)
## rectx - rectangles (list of min and max points) 
##
## Returns
## list of min and max points of the grid for each rectangle 
######################################################################################

find.gridpts <- function(gridx, suppx)
{
  xmax <- suppx$xmax
  xmin <- suppx$xmin
  d <- ncol(xmax)
  n <- nrow(xmax)
  gridpts.min <- matrix(0, ncol=d, nrow=n)
  gridpts.max <- gridpts.min
  
  for (i in 1:n)
    for (j in 1:d)    
    {
      ## find index of last element of gridx smaller than min support  
      tsum <- sum(xmin[i,j] >= gridx[[j]])
      if (tsum==0)
        gridpts.min[i,j] <- 1
      else
        gridpts.min[i,j] <- tsum

      ## find index of first element gridx greater than max support 
      gridpts.max[i,j] <- sum(xmax[i,j] >= gridx[[j]])
    }   
        
  return(list(xmin=gridpts.min, xmax=gridpts.max))
} 

######################################################################################
## Find the nearest grid points surrounding point x
######################################################################################

find.nearest.gridpts <- function(x, gridx, f)
{
  if (!is.list(gridx))
    return(find.nearest.gridpts.1d(x=x, gridx=gridx, f=f))
  else
  {  
    if (is.vector(x)) x <- as.matrix(t(x))
    d <- ncol(x)
    n <- nrow(x)
    gridsize <- sapply(gridx,length)
    gind <- matrix(0, nrow=n, ncol=d)
    
    for (i in 1:n)
      for (j in 1:d)
      {
        tsum <- sum(x[i,j] >= gridx[[j]])
        if (tsum==0)
          gind[i,j] <- 1
        else
          gind[i,j] <- tsum
      }
  }

  bperm <- list()
  for (j in 1:d) bperm[[j]] <- elem(1,2)
  binary.perm <- as.matrix(expand.grid(bperm))
  colnames(binary.perm) <- NULL

  gind.list <- list()
  fx <- rep(0, length=n)
  for (i in 1:n)
  {
    gind.list[[i]] <- matrix(gind[i,], nrow=2^d, ncol=d, byrow=TRUE) + binary.perm
    w <- matrix(0, nrow=2^d, ncol=d)
    gridw <- matrix(0, nrow=2^d, ncol=d)
    for (j in 1:d)
    {  
      gind.list[[i]][gind.list[[i]][,j]>=gridsize[j]] <- gridsize[j]
      gridw[,j] <- gridx[[j]][gind.list[[i]][,j]]
    }
    w <- 1/apply((matrix(x[i,], nrow=2^d, ncol=d, byrow=TRUE) - gridw)^2, 1, sum)
    w[w>1e5] <- 1e5
    w <- w/sum(w)
    fx[i] <- sum(w*f[gind.list[[i]]])
  }

  return(list(index=gind.list, fx=fx))
}


find.nearest.gridpts.1d <- function(x, gridx, f)
{
  n <- length(x)
  gind <- rep(0, length=n)

  for (i in 1:length(x))
  {
    tsum <- sum(x[i] >= gridx)
    if (tsum==0)
      gind[i] <- 1
    else
      gind[i] <- tsum
  }

  gind2 <- gind+1
  gind2[gind2>length(gridx)] <- length(gridx)
  gind <- cbind(gind, gind2)
  colnames(gind) <- NULL

  fx <- rep(0, n)
  for (i in 1:n)
  {
    w <- 1/(x[i] - gridx[gind[i,]])^2
    w[w>1e5] <- 1e5
    w <- w/sum(w)
   
    fx[i] <- sum(w*f[gind[i,]])
  }

  return(list(index=gind, fx=fx))
}


#########################################################################################################
## Multivariate kernel density estimate using normal kernels
##
## Parameters
## x - points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
## eval.points - compute density estimate at these points (if missing
##            and dim(x) = 2, 3 compute density estimate over grid)  
## eval.levels - compute 3-D in 2-D slices taken at these level curves   
##
## Returns
## list with first d components with the points that the density
## estimate is evaluated at, and values of the density estimate 
#############################################################################################################


kde <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, compute.cont=FALSE, approx.cont=TRUE)
{
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
    
    if (positive & d==1)
    {
      y <- log(x)
      fhat <- kdde.binned(x=y, H=H, h=h, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.order=0)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
    {
      ##if (d>1){ if (!identical(diag(diag(H)), H)) stop("Binned estimation defined for diagonal H only")}
      fhat <- kdde.binned(x=x, H=H, h=h, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.order=0)
    }
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
      {
        fhat <- kde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, xmin=xmin, xmax=xmax, adj.positive=adj.positive, gridtype=gridtype, w=w)
      }
      else
        fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive, w=w)
    }
     ## multi-dimensional
     else
     {  
       if (is.data.frame(x)) x <- as.matrix(x)

       if (missing(eval.points))
       {
         if (d==2)
           fhat <- kde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w)
         else if (d == 3)
           fhat <- kde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w) 
         else 
           stop("Need to specify eval.points for more than 3 dimensions")
       }
       else
         fhat <- kde.points(x, H, eval.points, w=w)     
     }
  }
  
  fhat$binned <- binned
  
  ## add variable names
  if (d==1)
  {
    x.names <- deparse(substitute(x))
  }
  else
  {  
    x.names <- colnames(x)
    if (is.null(x.names))
    {
      x.names <- strsplit(deparse(substitute(x)), "\\[")[[1]][1]
      x.names <- paste(x.names, "[,", 1:d,"]",sep="") 
    }
  }
  fhat$names <- x.names
  fhat$w <- w
  class(fhat) <- "kde"

  ## compute prob contour levels
  if (compute.cont & missing(eval.points))
    fhat$cont <- contourLevels(fhat, cont=1:99, approx.cont=approx.cont)

  return(fhat)
 }

#########################################################################################
## Univariate kernel density estimate on a grid
#########################################################################################

kde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype, w)
{
  if (missing(xmin)) xmin <- min(x) - h*supp
  if (missing(xmax)) xmax <- max(x) + h*supp
  if (missing(gridtype)) gridtype <- "linear"
  
  if (positive)
  {
    if (missing(adj.positive)) adj.positive <- abs(min(x))
    y <- log(x + adj.positive)  ## transform positive data x to real line

    gridx <- seq(max(0, xmin), xmax, length=gridsize)
    gridy <- log(gridx + adj.positive)
    gridtype.vec <- "linear" 
  }
  else
  {
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
  }
  n <- length(y)
  
  est <- dnorm.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=w/n)
  fhat <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE)
  
  if (positive)
  {
    ## compute transformation KDE
    fhat$estimate <- fhat$estimate / exp(gridy)
    fhat$x <- x
    fhat$eval.points <- gridx 
  }
  
  class(fhat) <- "kde"
  
  return(fhat)
}


#####################################################################################
## Bivariate kernel density estimate using normal kernels, evaluated over grid
##
## Parameters
## x - data points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
######################################################################################

kde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w)
{
  ## initialise grid 
  n <- nrow(x)
  if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 

  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
  
  for (i in 1:n)
  {
    ## compute evaluation points 
    eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
    eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
    eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
    eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
    eval.x.len <- length(eval.x)
    eval.pts <- permute(list(eval.x, eval.y))
    fhat <- dmvnorm(eval.pts, x[i,], H)
    
    ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
    for (j in 1:length(eval.y))
      fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
        fhat.grid[eval.x.ind, eval.y.ind[j]] + 
          w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
  }
  
  fhat.grid <- fhat.grid/n
  gridx1 <- list(gridx[[1]], gridx[[2]]) 
  
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE)
  
  return(fhat.list)
}


#######################################################################################
## Trivariate kernel density estimate using normal kernels, evaluated over grid
##
## Parameters
## x - data points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
#######################################################################################

kde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w)
{
  ## initialise grid 
  n <- nrow(x)

  if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), 
               length(gridx[[3]])))
  
  for (i in 1:n)
  {
    ## compute evaluation points
    eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
    eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
    eval.z <- gridx[[3]][grid.pts$xmin[i,3]:grid.pts$xmax[i,3]]
    eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
    eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
    eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
    eval.x.len <- length(eval.x)
    eval.pts <- permute(list(eval.x, eval.y))
   
    ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 

    for (k in 1:length(eval.z))
    {
      fhat <- w[i]*dmvnorm(cbind(eval.pts, eval.z[k]), x[i,], H)
      for (j in 1:length(eval.y))
        fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
          fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
            fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
     }
  }
  
  fhat.grid <- fhat.grid/n

  gridx1 <- list(gridx[[1]], gridx[[2]], gridx[[3]]) 
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE)

  return(fhat.list)
}


##########################################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
##
## Parameters
## x - data points
## H - bandwidth matrix
## eval.points - points where to evaluate density estimate
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
##########################################################################################

kde.points <- function(x, H, eval.points, w) 
{
  n <- nrow(x)
  Hs <- numeric(0)
  for (i in 1:n)
    Hs <- rbind(Hs, H)

  fhat <- dmvnorm.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n)

  return(list(x=x, eval.points=eval.points, estimate=fhat, H=H, gridded=FALSE))
}

kde.points.1d <- function(x, h, eval.points, positive=FALSE, adj.positive, w) 
{
  n <- length(x)

  if (positive)
  {
    if (missing(adj.positive)) adj.positive <- abs(min(x))
    y <- log(x + adj.positive)  ## transform positive data x to real line
    eval.pointsy <- log(eval.points + adj.positive)
  }
  else
  {
    y <- x
    eval.pointsy <- eval.points
  }
  
  fhat <- dnorm.mixt(x=eval.pointsy, mus=y, sigmas=rep(h,n), props=w/n)
  if (positive)
    fhat <- fhat/(eval.points + adj.positive) ##fhat/exp(eval.pointsy)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2, gridded=FALSE))
}

#### sum of KDE evaluated at eval.points

kde.points.sum <- function(x, H, eval.points, verbose=FALSE, binned=FALSE, bgridsize)
{
  nx <- nrow(x)
  ne <- nrow(eval.points)
  d <- ncol(x)
  
  if (binned)
  {
    fhatx <- kdde(x=x, deriv.order=0, H=H, binned=TRUE, bgridsize=bgridsize)
    fhat <- find.nearest.gridpts(x=eval.points, gridx=fhatx$eval.points, f=fhatx$estimate)$fx
    fhat.sum <- sum(fhat)
    fhat.sumsq <- sum(fhat^2)
  }
  else
  {
    if (verbose) pb <- txtProgressBar() 
    n.per.group <- max(c(round(1e6/sqrt(ne*nx)),1))
    ##ngroup <- max(ne%/%n.per.group+1,1)
    n.seq <- seq(1, ne, by=n.per.group)
    if (tail(n.seq,n=1) <= ne) n.seq <- c(n.seq, ne+1)

    fhat.sum <- 0
    fhat.sumsq <- 0
    if (length(n.seq)> 1)
    {
      for (i in 1:(length(n.seq)-1))
      {  
        if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1)) 
        difs <- differences(x=x, y=eval.points[n.seq[i]:(n.seq[i+1]-1),])
        fhat <- dmvnorm(x=difs, mean=rep(0,d), sigma=H)
        fhat <- apply(matrix(fhat, nrow=nx), 2, sum)/nx
        fhat.sum <- fhat.sum + sum(fhat)
        fhat.sumsq <- fhat.sumsq + sum(fhat^2)
      }
    }
    else
    {
      fhat <- dmvnorm(x=x, mean=eval.points, sigma=H)
      fhat <- apply(matrix(fhat, nrow=nx), 2, sum)/nx
      fhat.sum <- sum(fhat)
      fhat.sumsq <- sum(fhat^2)
    }
    if (verbose) close(pb)
  }
    
  return(list(sum=fhat.sum, sumsq=fhat.sumsq))
}


#######################################################################################
## Display kernel density estimate
##
## Parameters
## fhat - output from call to `kde'
#######################################################################################

plot.kde <- function(x, ...)
{ 
  fhat <- x

  if (is.vector(fhat$x))
    plotkde.1d(fhat, ...)
  else
  {
    d <- ncol(fhat$x)

    if (d==2) 
    {
      plotret <- plotkde.2d(fhat, ...)
      invisible(plotret)
    }
    else if (d==3)
    {
      plotkde.3d(fhat, ...)
      invisible()
    }
    else 
      stop ("Plot function only available for 1, 2 or 3-dimensional data")
  }
}

plotkde.1d <- function(fhat, xlab, ylab="Density function", add=FALSE,
  drawpoints=FALSE, ptcol="blue", jitter=FALSE, ...) 
{
  if (missing(xlab)) xlab <- fhat$names
  
  if (add)
    lines(fhat$eval.points, fhat$estimate, xlab=xlab, ylab=ylab, ...)
  else
    plot(fhat$eval.points, fhat$estimate, type="l", xlab=xlab, ylab=ylab, ...) 

  if (drawpoints)
    if (jitter)
      rug(jitter(fhat$x), col=ptcol)
    else
      rug(fhat$x, col=ptcol)
}


#########################################################################################
## Display bivariate kernel density estimate
##
## Parameters 
## fhat - output from 'kde.grid'
## display - "persp" - perspective plot
##         - "slice" - contour plot
##         - "image" image plot
## cont - vector of contours to be plotted
#########################################################################################

plotkde.2d <- function(fhat, display="slice", cont=c(25,50,75), abs.cont, approx.cont=FALSE,
    xlab, ylab, zlab="Density function", cex=1, pch=1, labcex,  
    add=FALSE, drawpoints=FALSE, drawlabels=TRUE, theta=-30, phi=40, d=4,
    ptcol="blue", col, lwd=1, ...) ##shade=0.75, border=NA, persp.col="grey", ...)
{
  disp1 <- substr(display,1,1)
  if (!is.list(fhat$eval.points))
    stop("Need a grid of density estimates")

  if (missing(xlab)) xlab <- fhat$names[1]
  if (missing(ylab)) ylab <- fhat$names[2]
  if (missing(labcex)) labcex <-1
  if (missing(approx.cont)) approx.cont <- (nrow(fhat$x) > 2000)

  ## perspective/wireframe plot
  if (disp1=="p")
  {
    hts <- seq(0, 1.1*max(fhat$estimate), length=100)
    if (missing(col)) col <- c("white", rev(heat.colors(length(hts))))
    if (length(col)<100) col <- rep(col, length=100)
    z <- fhat$estimate
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, length(hts)+1)
    plotret <- persp(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, col=col[facetcol], ...)
          ##shade=shade, border=border, col=persp.col, ...)
  }
  else if (disp1=="s") 
  {
    if (!add)
      plot(fhat$x[,1], fhat$x[,2], type="n",xlab=xlab, ylab=ylab, ...)

    ## compute contours
    if (missing(abs.cont))
    {
      if (!is.null(fhat$cont))
      {
        cont.ind <- rep(FALSE, length(fhat$cont))
        for (j in 1:length(cont))
          cont.ind[which(cont[j] == as.numeric(unlist(strsplit(names(fhat$cont),"%"))))] <- TRUE
        
        if (all(!cont.ind))
          hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        else
          hts <- fhat$cont[cont.ind]
      }
      else
        hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
    }
    ##else if (is.null(abs.cont))
    ##  hts <- contourLevels(fhat, n.pretty=5)  
    else
      hts <- abs.cont 
    
    hts <- sort(hts)
    
    if (missing(col)) col <- 1 #rev(heat.colors(length(hts)))
    if (length(col)<length(hts)) col <- rep(col, times=length(hts))
    
    ## draw contours         
    for (i in 1:length(hts)) 
    {
      if (missing(abs.cont)) scale <- cont[i]/hts[i]
      else scale <- 1

      if (hts[i]>0)
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate*scale, level=hts[i]*scale, add=TRUE, drawlabels=drawlabels, labcex=labcex, col=col[i], lwd=lwd, ...)
    }
 
    ## add points 
    if (drawpoints)
      points(fhat$x[,1], fhat$x[,2], col=ptcol, cex=cex, pch=pch)
  }
  ## image plot
  else if (disp1=="i")
  {
    image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
            xlab=xlab, ylab=ylab, add=add, ...)
    box()
  }
  else if (disp1=="f")
  {
    if (display=="filled.contour2")
    {
      ## compute contours
      if (missing(abs.cont))
      {
        if (!is.null(fhat$cont))
        {
          cont.ind <- rep(FALSE, length(fhat$cont))
          for (j in 1:length(cont))
            cont.ind[which(cont[j] == as.numeric(unlist(strsplit(names(fhat$cont),"%"))))] <- TRUE
          
          if (all(!cont.ind))
            hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
          else
            hts <- fhat$cont[cont.ind]
        }
        else
          hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
      }  
      else
        hts <- abs.cont 
      
      hts <- sort(hts)
      
      if (missing(col)) col <- c("transparent", rev(heat.colors(length(hts))))
      
      clev <- c(-0.01*max(abs(fhat$estimate)), hts, max(c(fhat$estimate, hts)) + 0.01*max(abs(fhat$estimate)))
      image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, xlab=xlab, ylab=ylab, add=add, col=col[1:(length(hts)+1)], breaks=clev, ...)

      ## draw contours         
     
      for (i in 1:length(hts)) 
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, level=hts[i], add=TRUE, drawlabels=FALSE, col=col[i+1], lwd=7)
      if (!missing(lwd))
      {
        for (i in 1:length(hts)) 
        {
          if (missing(abs.cont)) scale <- cont[i]/hts[i]
          else scale <- 1
          
          if (lwd >=1) contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate*scale, level=hts[i]*scale, add=TRUE, drawlabels=drawlabels, col=1, labcex=labcex, lwd=lwd, ...)
        }
      }
    }
    else
       filled.contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, xlab=xlab, ylab=ylab, ...)
  }
  if (disp1=="p")  invisible(plotret)
  else invisible()
}
  


#######################################################################################
## Display trivariate kernel density estimate
#######################################################################################


plotkde.3d <- function(fhat, cont=c(25,50,75), abs.cont, approx.cont=FALSE, colors, alphavec, size=3, ptcol="blue", add=FALSE, 
  xlab, ylab, zlab, drawpoints=FALSE, alpha=1, box=TRUE, axes=TRUE, ...)

{
  ##require(rgl)
  ##require(misc3d)
  
  if (missing(approx.cont))
    approx.cont <- (nrow(fhat$x) > 2000)
    
  ##if (missing(abs.cont))
  ##  hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
  ##else
  ##  hts <- abs.cont

  ## compute contours
  if (missing(abs.cont))
  {
    if (!is.null(fhat$cont))
      {
        cont.ind <- rep(FALSE, length(fhat$cont))
          for (j in 1:length(cont))
            cont.ind[which(cont[j] == as.numeric(unlist(strsplit(names(fhat$cont),"%"))))] <- TRUE
          
        if (all(!cont.ind))
          hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        else
          hts <- fhat$cont[cont.ind]
      }
    else
      hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
  }  
  else
    hts <- abs.cont
  
  nc <- length(hts)
  
  if (missing(colors)) colors <- rev(heat.colors(nc))
  if (missing(xlab)) xlab <- fhat$names[1]
  if (missing(ylab)) ylab <- fhat$names[2]
  if (missing(zlab)) zlab <- fhat$names[3]
  if (missing(alphavec)) alphavec <- seq(0.1,0.5,length=nc)

  if (drawpoints)
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=ptcol, alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
  else
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], type="n", xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
  bg3d(col="white")
  
  for (i in 1:nc)
    if (hts[nc-i+1] < max(fhat$estimate))
      contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE, color=colors[i], alpha=alphavec[i], box=FALSE, axes=FALSE, ...)

  if (box) box3d()
  if (axes) axes3d()
}




######################################################################################
## Contour levels 
######################################################################################

## create S3 generic 
contourLevels <- function(x, ...){  
  UseMethod("contourLevels")  
}   

contourLevels.kde <- function(x, prob, cont, nlevels=5, approx=FALSE, ...)
{ 
  fhat <- x
  if (is.vector(fhat$x))
  {
    d <- 1; n <- length(fhat$x)
  }
  else
  {
    d <- ncol(fhat$x); n <-nrow(fhat$x)
    if (!is.matrix(fhat$x)) fhat$x <- as.matrix(fhat$x)
  }

  if (is.null(x$w)) w <- rep(1, n)
  else w <- x$w

  if (is.null(fhat$gridded))
  {
    if (d==1) fhat$gridded <- fhat$binned
    else fhat$gridded <- is.list(fhat$eval.points)
  }
  
  if (missing(prob) & missing(cont))
    hts <- pretty(x$estimate, n=nlevels) 
  else
  {
    if (approx & fhat$gridded)
      dobs <- find.nearest.gridpts(x=fhat$x, gridx=fhat$eval.points, f=fhat$estimate)$fx
    else
      dobs <- kde(x=fhat$x, H=fhat$H, eval.points=fhat$x, w=w)$estimate 
    
    if (!missing(prob) & missing(cont))
      hts <- quantile(dobs, prob=prob)
    
    if (missing(prob) & !missing(cont))
      hts <- quantile(dobs, prob=(100-cont)/100)
  }
  
  return(hts)
}

