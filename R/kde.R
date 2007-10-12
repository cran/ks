
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

make.grid.ks <- function(x, H, tol, gridsize)
{
  d <- ncol(x)
  tol.H <-  tol * diag(H)
  minx <- apply(x, 2, min) - tol.H
  maxx <- apply(x, 2, max) + tol.H
  stepsize <- rep(0, d)
  gridx <- numeric(0)
 
  for (i in 1:d)
  {
    gridx <- c(gridx, list(seq(minx[i], maxx[i], length=gridsize[i])))
    stepsize[i] <- abs(gridx[[i]][1] - gridx[[i]][2])   
  }
  gridx <- c(gridx, list(stepsize = stepsize))
  
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
  minx <- matrix(0, nr=n, nc=d)
  maxx <- matrix(0, nr=n, nc=d)

  for (i in 1:n)
  {
    minx[i,] <- x[i,] - tol.H
    maxx[i,] <- x[i,] + tol.H 
  }
           
  return(list(minx = minx, maxx = maxx))
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
  maxx <- suppx$maxx
  minx <- suppx$minx
  d <- ncol(maxx)
  n <- nrow(maxx)
  gridpts.min <- matrix(0, nc=d, nr=n)
  gridpts.max <- gridpts.min
  
  for (i in 1:n)
    for (j in 1:d)    
    {
      # find index of last element of gridx smaller than min support  
      tsum <- sum(minx[i,j] >= gridx[[j]])
      if (tsum==0)
        gridpts.min[i,j] <- 1
      else
        gridpts.min[i,j] <- tsum

      # find index of first element gridx greater than max support 
      gridpts.max[i,j] <- sum(maxx[i,j] >= gridx[[j]])
    }   
        
  return(list(minx=gridpts.min, maxx=gridpts.max))
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


kde <- function(x, H, h, gridsize, binned=FALSE, bgridsize, supp=3.7, eval.points)
{ 
  if (binned)
  {
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
    fhat <- kde.binned(x=x, H=H, h=h, bgridsize=bgridsize, supp=supp)
  }
  else
  {
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
  
    
    if (is.vector(x))
    {
      if (!missing(h))
        fhat <- kde.1d(x=x, h=h, gridsize=gridsize, eval.points=eval.points)
      else if (!missing(H))
        fhat <- kde.1d(x=x, h=sqrt(H), gridsize=gridsize, eval.points=eval.points)
    }
    else
    {  
      if (is.data.frame(x)) x <- as.matrix(x)
                                                
      if (missing(eval.points))
      {
        d <- ncol(x)
        if (d==2)
          fhat <- kde.grid.2d(x, H, gridsize, supp)
        else if (d == 3)
          fhat <- kde.grid.3d(x, H, gridsize, supp) 
        else 
          stop("Need to specify eval.points for more than 3 dimensions")
      }
      else
        fhat <- kde.points(x, H, eval.points)
    }

  }

  fhat$binned <- binned
  class(fhat) <- "kde"
 
  
  return(fhat)
}

###############################################################################
### Multivariate binned kernel density estimate using normal kernels
###############################################################################


kde.binned <- function(x, H, h, bgridsize, supp)
{
  if (is.vector(x))
  {
    d <- 1
    n <- length(x)
  }
  else
  {
    d <- ncol(x)
    n <- nrow(x)
  }
  RK <- (4*pi)^(-d/2)

  if (d > 4)
    stop("Binning only available for 1- to 4-dim data")

  if (!identical(diag(diag(H)), H))
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
  bin.par <- dfltCounts.ks(x, bgridsize, sqrt(diag(H)), supp=3.7)  
  
  if (missing(h))
    fhat.grid <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=sqrt(diag(H)),
                           binned=TRUE, range.x=bin.par$range.x, se=FALSE)
  else 
    fhat.grid <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=h,
                           binned=TRUE, range.x=bin.par$range.x, se=FALSE)

  eval.points <- fhat.grid$x.grid
  fhat.grid <- fhat.grid$est
  fhat.grid[fhat.grid<0] <- 0
  #for (j in 1:length(bin.par$range.x))
  #  eval.points[[j]] <- seq(min(bin.par$range.x[[j]]), max(bin.par$range.x[[j]]), length=bgridsize[j])
  
  if (missing(h))
    fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, H=H)
  else
    fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, H=diag(h^2))
  
  return(fhat)
  
}

###############################################################################
## Univariate
###############################################################################

kde.1d <- function(x, h, gridsize, eval.points)
{
 
  if (missing(eval.points))
  {
    fhat <- bkde(x=x, bandwidth=h, gridsize=gridsize)
    fhat <- list(x=x, eval.points=fhat$x, estimate=fhat$y, h=h)
  }
  else
    fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points)

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

kde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL)
{
  # initialise grid 
  n <- nrow(x)
  if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize) 
  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))

  for (i in 1:n)
  {
    ## compute evaluation points 
    eval.x <- seq(gridx[[1]][grid.pts$minx[i,1]], 
                  gridx[[1]][grid.pts$maxx[i,1]], by=gridx$stepsize[1])
    eval.y <- seq(gridx[[2]][grid.pts$minx[i,2]], 
                  gridx[[2]][grid.pts$maxx[i,2]], by=gridx$stepsize[2])
    eval.x.ind <- c(grid.pts$minx[i,1]:grid.pts$maxx[i,1])
    eval.y.ind <- c(grid.pts$minx[i,2]:grid.pts$maxx[i,2])
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
  
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H)
  
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


kde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL)
{
  ## initialise grid 
  n <- nrow(x)

  if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize) 
  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), 
               length(gridx[[3]])))
  
  for (i in 1:n)
  {
    ## compute evaluation points 
    eval.x <- seq(gridx[[1]][grid.pts$minx[i,1]], 
                  gridx[[1]][grid.pts$maxx[i,1]], by=gridx$stepsize[1])
    eval.y <- seq(gridx[[2]][grid.pts$minx[i,2]], 
                  gridx[[2]][grid.pts$maxx[i,2]], by=gridx$stepsize[2])
    eval.z <- seq(gridx[[3]][grid.pts$minx[i,3]], 
                  gridx[[3]][grid.pts$maxx[i,3]], by=gridx$stepsize[3])
 
    eval.x.ind <- c(grid.pts$minx[i,1]:grid.pts$maxx[i,1])
    eval.y.ind <- c(grid.pts$minx[i,2]:grid.pts$maxx[i,2])
    eval.z.ind <- c(grid.pts$minx[i,3]:grid.pts$maxx[i,3])
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
  fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H)

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

kde.points.1d <- function(x, h, eval.points) 
{
    n <- length(x)
    fhat <- dnorm.mixt(x=eval.points, mus=x, sigmas=rep(h, n), props=rep(1, n)/n)

    return(list(x=x, eval.points=eval.points, estimate=fhat, h=h))
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
    plotkde.1d(fhat, drawpoints=drawpoints, ...)
  else
  {
    d <- ncol(fhat$x)

    if (d==2) 
      plotkde.2d.v2(fhat, drawpoints=drawpoints, ...)
    else if (d==3)
      plotkde.3d(fhat, drawpoints=drawpoints, ...)
    else 
      stop ("Plot function only available for 1, 2 or 3-dimensional data")
  }
}

plotkde.1d <- function(fhat, xlab="x", ylab="Density function", add=FALSE,
  drawpoints=TRUE, ptcol="blue", lcol="black", ...)
{
  if (add)
    lines(fhat$eval.points, fhat$estimate, col=lcol, xlab=xlab, ylab=ylab, ...)
  else
    plot(fhat$eval.points, fhat$estimate, col=lcol, type="l", xlab=xlab, ylab=ylab, ...)

  if (drawpoints)
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

plotkde.2d.v2 <- function(fhat, display="slice", cont=c(25,50,75), abs.cont, cex=0.7, 
    xlab, ylab, zlab="Density function",  
    add=FALSE, drawpoints=TRUE, drawlabels=TRUE, theta=-30, phi=40, d=4,
    pch, ptcol="blue", lcol="black",
    ...)
{
  disp1 <- substr(display,1,1)
  if (!is.list(fhat$eval.points))
    stop("Need a grid of density estimates")

  if (is.data.frame(fhat$x))
  {
    xlab <- names(fhat$x)[1]
    ylab <- names(fhat$x)[2]
  }

  x.names <- colnames(fhat$x) 
  if (!is.null(x.names))
  {
    if (missing(xlab)) xlab <- x.names[1]
    if (missing(ylab)) ylab <- x.names[2]
  }
  else
  {
    xlab <- "x"
    ylab <- "y"
  }
 
  eval1 <- fhat$eval.points[[1]]
  eval2 <- fhat$eval.points[[2]]
  
  ## perspective/wireame plot
  if (disp1=="p")
    persp(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
          theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, ...)
  
  else if (disp1=="s")
  {
    if (!add)
      plot(fhat$x[,1], fhat$x[,2], type="n",xlab=xlab, ylab=ylab, ...)

    n <- nrow(fhat$x)
    RK <- (4*pi)^(-2/2)
    bgridsize <- dim(fhat$estimate)
    
    if (missing(abs.cont))
    {
      ## for large sample sizes, use binned approx. 
      if (fhat$binned | nrow(fhat$x) >= 5e3)
      {
        bin.par <- dfltCounts.ks(fhat$x, bgridsize, sqrt(diag(fhat$H)), supp=3.7)
        dobs <- rep(fhat$estimate, round(bin.par$counts,0))
      }
      else
        dobs <- kde(x=fhat$x, H=fhat$H, eval.points=fhat$x)$estimate 
      
      hts <- quantile(dobs, prob=(100 - cont)/100)
    }
    else
      hts <- abs.cont
    
    ## compute and draw contours

    if (missing(abs.cont))
      for (i in 1:length(hts)) 
      {
        scale <- cont[i]/hts[i]
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                fhat$estimate*scale, level=hts[i]*scale, add=TRUE, 
                drawlabels=drawlabels, col=lcol, ...)
      }
    else
      contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
              fhat$estimate, level=hts, add=TRUE, 
              drawlabels=drawlabels, col=lcol, ...)
    

    ## add points 
    if (drawpoints)
      points(fhat$x[,1], fhat$x[,2], cex=cex, col=ptcol)
  }
  ## image plot
  else if (disp1=="i")
  {
    image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
            xlab=xlab, ylab=ylab, ...)
    box()
  }
  else if (disp1=="f")
    filled.contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
                   xlab=xlab, ylab=ylab, ...)
}
  




###############################################################################
## Display trivariate kernel density estimate
###############################################################################


plotkde.3d <- function(fhat, cont=c(25,50,75), abs.cont, colors,
  alphavec, size=3, ptcol="blue", add=FALSE, 
  xlab, ylab, zlab, drawpoints=FALSE, ...)

{
  n <- nrow(fhat$x)
  RK <- (4*pi)^(-3/2)
  bgridsize <- dim(fhat$estimate)
  
  ## for large sample sizes, use binned approx. 
  if (missing(abs.cont))
  {
    if (fhat$binned | nrow(fhat$x) >= 5e4)
    {
      bin.par <- dfltCounts.ks(fhat$x, bgridsize, sqrt(diag(fhat$H)), supp=3.7)
      dobs <- rep(fhat$estimate, round(bin.par$counts,0))
    }
    else
      dobs <- kde(x=fhat$x, H=fhat$H, eval.points=fhat$x)$estimate 
  
    hts <- quantile(dobs, prob = (100-cont)/100)
  }
  else
    hts <- abs.cont
  
  nc <- length(hts)
  
  if (missing(colors))
    colors <- rev(heat.colors(nc))

 
  x.names <- colnames(fhat$x)

  if (missing(xlab))
    if (is.null(x.names)) xlab <- "x" else xlab <- x.names[1]
  if (missing(ylab))
    if (is.null(x.names)) ylab <- "y" else ylab <- x.names[2]
  if (missing(zlab))
    if (is.null(x.names)) zlab <- "z" else zlab <- x.names[3]
  
  if (missing(alphavec))
    alphavec <- seq(0.1,0.5,length=nc)

  
  if (!add) clear3d()
  
  bg3d(col="white")
  if (drawpoints)
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=ptcol, alpha=1, xlab=xlab, ylab=ylab, zlab=zlab, ...)
  else
    plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], type="n", xlab=xlab, ylab=ylab, zlab=zlab, ...)
  
  for (i in 1:nc)
  {
    ##scale <- cont[i]/hts[i]
    contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]],
              y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE,
              color=colors[i], alpha=alphavec[i], ...)
  }
}


### highest density regions

prob.cont <- function(x, H, prob=c(0.25,0.50,0.75), binned=FALSE, bgridsize)
{
  if (is.vector(x)) d <-1
  else d <- ncol(x)
  
  if (missing(bgridsize))
    if (d==1)
      bgridsize <- 401
    else if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)
  
  if (binned)
  {
    fhat <- kde(x=x, H=H, binned=binned)
    bin.par <- dfltCounts.ks(x, bgridsize, sqrt(diag(fhat$H)), supp=3.7)
    fhat.est <- rep(fhat$estimate, round(bin.par$counts,0))
  }
  else
    fhat.est <- kde.points(x=x, H=H, eval.points=x)$est
  
  cont <- quantile(fhat.est, prob=prob)

  return(cont)
}
