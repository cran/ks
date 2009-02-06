###############################################################################
# Multivariate kernel density estimators
###############################################################################


###############################################################################
# Generate grid over a set of points
#
# Parameters
# x - data points
# H - bandwidth matrix
# tol - tolerance = extra coverage exceeding the range of x   
# gridsize - number of points for each direction
#
# Returns
# gridx - list of intervals, one for each co-ord direction so that
#         gridx[[1]] x gridx[[2]] x ... x gridx[[d]] is the grid
# stepsize - vector of step sizes 
###############################################################################

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


###############################################################################
# Generate kernel (rectangular) support at data point
# 
# Parameters
# x - data points
# H - bandwidth matrix
# tol - tolerance = extra coverage exceeding the range of x 
#
# Returns
# list of min and max points of support (here we parameterise rectangles
# by their min = lower left co-ord and max = upper right coord)
###############################################################################

make.supp <- function(x, H, tol)
{
  n <- nrow(x)
  d <- ncol(x)
  tol.H <- tol * diag(H)
  xmin <- matrix(0, nr=n, nc=d)
  xmax <- matrix(0, nr=n, nc=d)

  for (i in 1:n)
  {
    xmin[i,] <- x[i,] - tol.H
    xmax[i,] <- x[i,] + tol.H 
  }
           
  return(list(xmin = xmin, xmax = xmax))
}


###############################################################################
# Find the grid points contained in kernel support rectangles 
#
# Parameters
# gridx - grid (list of subdivided intervals)
# rectx - rectangles (list of min and max points) 
#
# Returns
# list of min and max points of the grid for each rectangle 
###############################################################################

find.gridpts <- function(gridx, suppx)
{
  xmax <- suppx$xmax
  xmin <- suppx$xmin
  d <- ncol(xmax)
  n <- nrow(xmax)
  gridpts.min <- matrix(0, nc=d, nr=n)
  gridpts.max <- gridpts.min
  
  for (i in 1:n)
    for (j in 1:d)    
    {
      # find index of last element of gridx smaller than min support  
      tsum <- sum(xmin[i,j] >= gridx[[j]])
      if (tsum==0)
        gridpts.min[i,j] <- 1
      else
        gridpts.min[i,j] <- tsum

      # find index of first element gridx greater than max support 
      gridpts.max[i,j] <- sum(xmax[i,j] >= gridx[[j]])
    }   
        
  return(list(xmin=gridpts.min, xmax=gridpts.max))
} 

###############################################################################
# Multivariate kernel density estimate using normal kernels
#
# Parameters
# x - points
# H - bandwidth matrix
# gridsize - number of interval points in grid
# supp - effective support of kernel
# eval.points - compute density estimate at these points (if missing
#            and dim(x) = 2, 3 compute density estimate over grid)  
# eval.levels - compute 3-D in 2-D slices taken at these level curves   
#
# Returns
# list with first d components with the points that the density
# estimate is evaluated at, and values of the density estimate 
###############################################################################


kde <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize,  positive=FALSE, adj.positive)
{
  ## compute binned estimator
  if (binned)
  {
    if (!missing(eval.points))
      stop("Both binned=TRUE and eval.points are non-empty.")
    
    if (missing(bgridsize))
    {
      if (is.vector(x)) 
        bgridsize <- 401
      else
      {
        d <- ncol(x)
        
        if (d==2)
          bgridsize <- rep(151,d)
        else if (d==3)
          bgridsize <- rep(51, d)
        else if (d==4)
          bgridsize <- rep(21, d)
      }
    }
    if (positive & is.vector(x))
    {
      y <- log(x)
      fhat <- kde.binned(x=y, H=H, h=h, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax)
      fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
      fhat$eval.points <- exp(fhat$eval.points)
      fhat$x <- x
    }
    else
      fhat <- kde.binned(x=x, H=H, h=h, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax)
  }
  else
  {
    ## compute exact (non-binned) estimator
    if (missing(gridsize))
    {  
      if (is.vector(x)) 
        gridsize <- 401
      else
      {
        d <- ncol(x)
        if (d==2)
          gridsize <- rep(151,d)
        else 
          gridsize <- rep(51, d)
      }
    }
  
    ## 1-dimensional    
    if (is.vector(x))
    {
      if (!missing(H) & !missing(h))
        stop("Both H and h are both specified.")

      if (missing(h))
        h <- sqrt(H)

      if (missing(eval.points))
        fhat <- kde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, xmin=xmin, xmax=xmax, adj.positive=adj.positive, gridtype=gridtype)
       else
         fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive)
     }
     ## multi-dimensional
     else
     {  
       if (is.data.frame(x)) x <- as.matrix(x)

       if (missing(eval.points))
       {
         d <- ncol(x)
         if (d==2)
           fhat <- kde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype)
         else if (d == 3)
           fhat <- kde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype) 
         else 
           stop("Need to specify eval.points for more than 3 dimensions")
       }
       else
         fhat <- kde.points(x, H, eval.points)
     }

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
 ### Multivariate binned kernel density estimate using normal kernels
 ###############################################################################


 kde.binned <- function(x, H, h, bgridsize, supp, xmin, xmax)
 {
   if (is.vector(x))
   {
     d <- 1
     ##n <- length(x)
   }
   else
   {
     d <- ncol(x)
     ##n <- nrow(x)
   }
   ##RK <- (4*pi)^(-d/2)

   if (d==1)
     if (missing(H)) H <- as.matrix(h^2)
     else {h <- sqrt(H); H <- as.matrix(H)}

   if (d > 4)
     stop("Binning only available for 1- to 4-dim data")

   if (!identical(diag(diag(H)), H) & d > 1)
     stop("Binning requires diagonal bandwidth matrix")

    
   if (missing(bgridsize))
     if (d==1)
       bgridsize <- 401
     else if (d==2)
       bgridsize <- rep(151,d)
     else if (d==3)
       bgridsize <- rep(51, d)
     else if (d==4)
       bgridsize <- rep(21, d)

   ## linear binning
   if (missing(xmin) | missing(xmax))
     bin.par <- dfltCounts.ks(x, bgridsize, sqrt(diag(H)), supp=3.7)  
   else
   { 
     xrange.list <- list()
     for (j in 1:d)
        xrange.list[[j]] <- c(xmin[j], xmax[j])

     bin.par <- dfltCounts.ks(x, bgridsize, sqrt(diag(H)), supp=3.7, range.x=xrange.list)  
   }

   fhat.grid <- drvkde(x=bin.par$counts, drv=rep(0,d),bandwidth=sqrt(diag(H)), binned=TRUE, range.x=bin.par$range.x, se=FALSE, gridsize=bgridsize)
   eval.points <- fhat.grid$x.grid
   fhat.grid <- fhat.grid$est
   fhat.grid[fhat.grid<0] <- 0

   if (d==1)
     fhat <- list(x=x, eval.points=unlist(eval.points), estimate=fhat.grid, H=h^2, h=h)
   else
     fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, H=H)

   return(fhat)

 }

 ###############################################################################
 ## Univariate kernel density estimate on a grid
 ###############################################################################

 kde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype)
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
       ##gridy <- seq(sign(xmin)*sqrt(abs(xmin)), sign(xmax)*sqrt(abs(xmax)), length=gridsize)^2
       gridtype.vec <- "sqrt"
     }
   }
   n <- length(y)
 
   est <- dnorm.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=rep(1,n)/n)
   fhat <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec)
 
   if (positive)
   {
     ## compute transformation KDE
     fhat$estimate <- fhat$estimate / exp(gridy)
     fhat$x <- x
     fhat$eval.points <- gridx # exp(gridy) - adj.positive
   }
   
   class(fhat) <- "kde"
   
   return(fhat)
}

###############################################################################
# Bivariate kernel density estimate using normal kernels, evaluated over grid
#
# Parameters
# x - data points
# H - bandwidth matrix
# gridsize - number of interval points in grid
# supp - effective support of kernel
#
# Returns
# list with fields
# x - data points
# eval.points - points that KDE is evaluated at
# estimate - KDE evaluated at eval.points 
# H - bandwidth matrix 
###############################################################################

kde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype)
{
  # initialise grid 
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
    ##eval.x <- seq(gridx[[1]][grid.pts$xmin[i,1]], 
    ##              gridx[[1]][grid.pts$xmax[i,1]], by=gridx$stepsize[1])
    ##eval.y <- seq(gridx[[2]][grid.pts$xmin[i,2]], 
    ##              gridx[[2]][grid.pts$xmax[i,2]], by=gridx$stepsize[2])

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
          fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
  }
  
  fhat.grid <- fhat.grid/n
  gridx1 <- list(gridx[[1]], gridx[[2]]) 
  
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype)
  
  return(fhat.list)
}


###############################################################################
# Trivariate kernel density estimate using normal kernels, evaluated over grid
#
# Parameters
# x - data points
# H - bandwidth matrix
# gridsize - number of interval points in grid
# supp - effective support of kernel
#
# Returns
# list with fields
# x - data points
# eval.points - points that KDE is evaluated at
# estimate - KDE evaluated at eval.points 
# H - bandwidth matrix 
###############################################################################


kde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype)
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
    ##eval.x <- seq(gridx[[1]][grid.pts$xmin[i,1]], 
    ##              gridx[[1]][grid.pts$xmax[i,1]], by=gridx$stepsize[1])
    ##eval.y <- seq(gridx[[2]][grid.pts$xmin[i,2]], 
    ##              gridx[[2]][grid.pts$xmax[i,2]], by=gridx$stepsize[2])
    ##eval.z <- seq(gridx[[3]][grid.pts$xmin[i,3]], 
    ##              gridx[[3]][grid.pts$xmax[i,3]], by=gridx$stepsize[3])
 
    eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
    eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
    eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
    eval.x.len <- length(eval.x)
    eval.pts <- permute(list(eval.x, eval.y))
   
    ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 

    for (k in 1:length(eval.z))
    {
      fhat <- dmvnorm(cbind(eval.pts, eval.z[k]), x[i,], H)
      for (j in 1:length(eval.y))
        fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
          fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
            fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
     }
  }
  
  fhat.grid <- fhat.grid/n

  gridx1 <- list(gridx[[1]], gridx[[2]], gridx[[3]]) 
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype)

  return(fhat.list)
}





###############################################################################
# Multivariate kernel density estimate using normal kernels,
# evaluated at each sample point
#
# Parameters
# x - data points
# H - bandwidth matrix
# eval.points - points where to evaluate density estimate
#
# Returns
# list with fields
# x - data points
# eval.points - points that KDE is evaluated at
# estimate - KDE evaluated at eval.points 
# H - bandwidth matrix 
###############################################################################

kde.points <- function(x, H, eval.points) 
{
  n <- nrow(x)
  Hs <- numeric(0)
  for (i in 1:n)
    Hs <- rbind(Hs, H)
  
  fhat <- dmvnorm.mixt(x=eval.points, mus=x, Sigmas=Hs, props=rep(1, n)/n)

  return(list(x=x, eval.points=eval.points, estimate=fhat, H=H))
}

kde.points.1d <- function(x, h, eval.points, positive=FALSE, adj.positive) 
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
  
  #fhat <- dnorm.mixt(x=eval.pointsy, mus=y, sigmas=rep(h, n), props=rep(1, n)/n)
  fhat <- 0
  for (k in 1:n)
    fhat <- fhat + dnorm(x=eval.pointsy, mean=y[k], sd=h)
  fhat <- fhat/n
  if (positive)
    fhat <- fhat/(eval.points + adj.positive) ##fhat/exp(eval.pointsy)
  
  return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2))
}



###############################################################################
# Display kernel density estimate
#
# Parameters
# fhat - output from call to `kde'
###############################################################################

plot.kde <- function(x, drawpoints=FALSE, ...)
{ 
  fhat <- x

  if (is.vector(fhat$x))
    plotkde.1d.v2(fhat, drawpoints=drawpoints, ...)
  else
  {
    d <- ncol(fhat$x)

    if (d==2) 
    {
      plotret <- plotkde.2d.v2(fhat, drawpoints=drawpoints, ...)
      invisible(plotret)
    }
    else if (d==3)
    {
      plotkde.3d(fhat, drawpoints=drawpoints, ...)
       invisible()
    }
    else 
      stop ("Plot function only available for 1, 2 or 3-dimensional data")
  }
}

plotkde.1d.v2 <- function(fhat, xlab, ylab="Density function", add=FALSE,
  drawpoints=TRUE, ptcol="blue", jitter=FALSE, ...) #col="black", ...)
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


###############################################################################
# Display bivariate kernel density estimate
#
# Parameters 
# fhat - output from 'kde.grid'
# display - "persp" - perspective plot
#         - "slice" - contour plot
#         - "image" image plot
# cont - vector of contours to be plotted
###############################################################################

plotkde.2d.v2 <- function(fhat, display="slice", cont=c(25,50,75), abs.cont,
    xlab, ylab, zlab="Density function", cex=1, pch=1,  
    add=FALSE, drawpoints=TRUE, drawlabels=TRUE, theta=-30, phi=40, d=4,
    ptcol="blue", ...) #shade=0.75, border=NA, persp.col="grey", ...)
{
  disp1 <- substr(display,1,1)
  if (!is.list(fhat$eval.points))
    stop("Need a grid of density estimates")

  if (missing(xlab)) xlab <- fhat$names[1]
  if (missing(ylab)) ylab <- fhat$names[2]

  ##eval1 <- fhat$eval.points[[1]]
  ##eval2 <- fhat$eval.points[[2]]
  
  ## perspective/wireframe plot
  if (disp1=="p")
    plotret <- persp(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
          theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, ...)
          ##shade=shade, border=border, col=persp.col, ...)
  
  else if (disp1=="s")
  {
    if (!add)
      plot(fhat$x[,1], fhat$x[,2], type="n",xlab=xlab, ylab=ylab, ...)
      
    ## compute contours
    if (missing(abs.cont))
      hts <- contourLevels(fhat, prob=(100-cont)/100)
    else if (is.null(abs.cont))
      hts <- contourLevels(fhat, n.pretty=5)  
    else
      hts <- abs.cont 
   
    ## draw contours         
    for (i in 1:length(hts)) 
    {
      if (missing(abs.cont)) scale <- cont[i]/hts[i]
      else scale <- 1

      contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
              fhat$estimate*scale, level=hts[i]*scale, add=TRUE, 
              drawlabels=drawlabels, ...)
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
    filled.contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
                   xlab=xlab, ylab=ylab, ...)

  if (disp1=="p")  invisible(plotret)
  else invisible()
}
  


###############################################################################
## Display trivariate kernel density estimate
###############################################################################


plotkde.3d <- function(fhat, cont=c(25,50,75), abs.cont, colors,
                       alphavec, size=3, ptcol="blue", add=FALSE, 
                       xlab, ylab, zlab, drawpoints=FALSE, alpha=1, ...)

{
  require(rgl)
  require(misc3d)
 
  if (missing(abs.cont))
    hts <- contourLevels(fhat, prob=(100-cont)/100)
  else
    hts <- abs.cont 
  nc <- length(hts)
  
  if (missing(colors))
    colors <- rev(heat.colors(nc))

  if (missing(xlab)) xlab <- fhat$names[1]
  if (missing(ylab)) ylab <- fhat$names[2]
  if (missing(zlab)) zlab <- fhat$names[3]
    #if (is.null(x.names)) zlab <- "z" else zlab <- x.names[3]
  
  if (missing(alphavec))
    alphavec <- seq(0.1,0.5,length=nc)

  
  ##if (!add) clear3d()
  
  bg3d(col="white")
  
  if (drawpoints)
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=ptcol, alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add, ...)
  else
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], type="n", xlab=xlab, ylab=ylab, zlab=zlab, add=add, ...)
  
  for (i in 1:nc)
    contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]],
              y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE,
              color=colors[i], alpha=alphavec[i], ...)
}




###############################################################################
### Contour levels 
###############################################################################

## create S3 generic 
contourLevels <- function(x, ...){  
  UseMethod("contourLevels")  
}   

contourLevels.kde <- function(x, prob, cont, nlevels=5, ...)
{
  fhat <- x
  if (is.vector(fhat$x))
  {
    d <- 1; n <- length(fhat$x); H <- as.matrix(fhat$H)
    bgridsize <- length(fhat$estimate)
  }
  else
  {
    d <- ncol(fhat$x); n <-nrow(fhat$x); H <- fhat$H
    bgridsize <- dim(fhat$estimate)
  }

  ## for large sample sizes, use binned approx. 
  if (n >= 5e3 & d <= 4 & fhat$binned)
  {
    bin.par <- dfltCounts.ks(fhat$x, bgridsize, sqrt(diag(H)), supp=3.7)
    dobs <- rep(fhat$estimate, round(bin.par$counts,0))
    dobs <- dobs[dobs>0]
  }
  else
    dobs <- kde(x=fhat$x, H=fhat$H, eval.points=fhat$x)$estimate 
  
  if (missing(prob) & missing(cont))
    hts <- pretty(dobs, n=nlevels) 
  
  if (!missing(prob) & missing(cont))
    hts <- quantile(dobs, prob=prob)
  
  if (missing(prob) & !missing(cont))
    hts <- quantile(dobs, prob=(100-cont)/100)
  
  return(hts)
}


