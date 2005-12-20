
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

make.grid <- function(x, H, tol, gridsize)
{
  d <- ncol(x)
  tol.H <-  tol * (diag(H) + abs(H[1,2])) 
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
  tol.H <- tol * (diag(H) + abs(H[1,2]))
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


kde <- function(x, H, gridsize, supp=3.7, eval.points)
{ 
  d <- ncol(x)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.list(x))
  {
    if (missing(gridsize))
      gridsize <- rep(50, d)

    if (missing(eval.points))
      fhat <- kde.pc.grid.2d(x, H, gridsize, supp)
    else
      fhat <- kde.pc.points(x, H, eval.points)
  }
  else
  {
    if (missing(gridsize))
      gridsize <- rep(50, d)

    if (missing(eval.points))
    {
      if (d == 2)
        fhat <- kde.grid.2d(x, H, gridsize, supp)
      else if (d == 3)
        fhat <- kde.grid.3d(x, H, gridsize, supp) 
      else 
        stop("Need to specify eval.points for more than 3 dimensions")
    }
    else
      fhat <- kde.points(x, H, eval.points)
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

kde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL)
{
  # initialise grid 
  n <- nrow(x)
  if (is.null(gridx))
    gridx <- make.grid(x, matrix.sqrt(H), tol=supp, gridsize=gridsize) 
  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))

  for (i in 1:n)
  {
    # compute evaluation points 
    eval.x <- seq(gridx[[1]][grid.pts$minx[i,1]], 
                  gridx[[1]][grid.pts$maxx[i,1]], by=gridx$stepsize[1])
    eval.y <- seq(gridx[[2]][grid.pts$minx[i,2]], 
                  gridx[[2]][grid.pts$maxx[i,2]], by=gridx$stepsize[2])
    eval.x.ind <- c(grid.pts$minx[i,1]:grid.pts$maxx[i,1])
    eval.y.ind <- c(grid.pts$minx[i,2]:grid.pts$maxx[i,2])
    eval.x.len <- length(eval.x)
    eval.pts <- permute(list(eval.x, eval.y))
    fhat <- dmvnorm(eval.pts, x[i,], H)
    
    # place vector of density estimate values `fhat' onto grid 'fhat.grid' 
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
  # initialise grid 
  n <- nrow(x)

  if (is.null(gridx))
    gridx <- make.grid(x, matrix.sqrt(H), tol=supp, gridsize=gridsize) 
  suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

  if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
  fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), 
               length(gridx[[3]])))
  
  for (i in 1:n)
  {
    # compute evaluation points 
    eval.x <- seq(gridx[[1]][grid.pts$minx[i,1]], 
                  gridx[[1]][grid.pts$maxx[i,1]], by=gridx$stepsize[1])
    eval.y <- seq(gridx[[2]][grid.pts$minx[i,2]], 
                  gridx[[2]][grid.pts$maxx[i,2]], by=gridx$stepsize[2])
    eval.z <- seq(gridx[[3]][grid.pts$minx[i,3]], 
                  gridx[[3]][grid.pts$maxx[i,3]], by=gridx$stepsize[3])
    #else
    #  eval.z <- eval.levels
 
    eval.x.ind <- c(grid.pts$minx[i,1]:grid.pts$maxx[i,1])
    eval.y.ind <- c(grid.pts$minx[i,2]:grid.pts$maxx[i,2])
    eval.z.ind <- c(grid.pts$minx[i,3]:grid.pts$maxx[i,3])
    eval.x.len <- length(eval.x)
    eval.pts <- permute(list(eval.x, eval.y))
   
    # place vector of density estimate values `fhat' onto grid 'fhat.grid' 

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
# Pre-clustered bivariate kernel density estimate using normal kernels
#
# Parameters
# x.pc - pre-clustered data points 
# H.pc - bandwidth matrices, stacked into a matrix
# gridsize - number of interval points in grid
# supp - effective support of kernel 
#
# Returns
# list with fields
# x - pre-clustered data points
# eval.points - points that KDE is evaluated at
# estimate - KDE evaluated at eval.points 
# H - bandwidth matrices
###############################################################################

kde.pc.grid.2d <- function(x.pc, H.pc, gridsize, supp=3.7)
{
  d <- ncol(x.pc$x)
  nu <- length(x.pc$nclust)
  x <- x.pc$x
  clust.ind <- x.pc$ind
  nclust <- x.pc$nclust

  # find largest bandwidth matrix to initialise grid
  detH <- vector() 
  for (j in 1:nu)
    detH[j] <- det(H.pc[((j-1)*d+1) : (j*d),])  
  Hmax.ind <- which.max(detH)
  Hmax <- H.pc[((Hmax.ind-1)*d+1) : (Hmax.ind*d),]

  # initialise grid 
  gridx <- make.grid(x, matrix.sqrt(Hmax), tol=supp, gridsize=gridsize) 
  fhat.grid <- array(0, dim = gridsize)
  n <- sum(nclust)
  
  for (j in 1:nu)
  {       
    xj <- x[clust.ind==j,]  # points in cluster j 
    if (is.vector(xj))
      xj <- t(matrix(xj))
    nj <- nclust[j]
    Hj <- H.pc[((j-1)*d+1) : (j*d),] 
    rectxj <- make.supp(xj, matrix.sqrt(Hj), tol=supp)
    grid.ptsj <- find.gridpts(gridx, rectxj)

    for (i in 1:nj)
    {
      # compute evaluation points
      eval.xj <- seq(gridx[[1]][grid.ptsj$minx[i,1]], 
                     gridx[[1]][grid.ptsj$maxx[i,1]], by=gridx$stepsize[1])
      eval.yj <- seq(gridx[[2]][grid.ptsj$minx[i,2]], 
                     gridx[[2]][grid.ptsj$maxx[i,2]], by=gridx$stepsize[2])
      eval.xj.ind <- c(grid.ptsj$minx[i,1]:grid.ptsj$maxx[i,1])
      eval.yj.ind <- c(grid.ptsj$minx[i,2]:grid.ptsj$maxx[i,2])
      eval.xj.len <- length(eval.xj)
      eval.pts <- permute(list(eval.xj, eval.yj))
      fhat <- dmvnorm(eval.pts, xj[i,], Hj) 

      # place vector of density estimate values `fhat' onto grid 'fhat.grid'
      for (k in 1:length(eval.yj))
        fhat.grid[eval.xj.ind, eval.yj.ind[k]] <- 
          fhat.grid[eval.xj.ind, eval.yj.ind[k]] + 
            fhat[((k-1) * eval.xj.len + 1):(k * eval.xj.len)]      
    }
  }

  fhat.grid <- fhat.grid/n
  gridx1 <- list(gridx[[1]], gridx[[2]])
    
  return(list(x=x.pc, eval.points=gridx1, estimate=fhat.grid, H=H.pc))
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


###############################################################################
# Preclustered multivariate kernel density estimate using normal kernels,
# evaluated at each sample point
#
# Parameters
# x.pc - pre-clustered data points
# Hs - bandwidth matrices
# eval.points - points where to evaluate density estimate
#
# Returns
# list with fields
# x - pre-clustered data points
# eval.points - points that KDE is evaluated at
# estimate - KDE evaluated at eval.points 
# H - bandwidth matrices  
###############################################################################

kde.pc.points <- function(x.pc, Hs, eval.points) 
{
  fhat <- 0
  d <- ncol(x.pc$x)
  nu <- length(x.pc$nclust)
  x <- x.pc$x
  clust.ind <- x.pc$ind
  nclust <- x.pc$nclust
  n <- sum(nclust)
  
  for (j in 1:nu)
  {  
    xj <- x[clust.ind==j,] # points in cluster j
    nj <- nclust[j]
    Hj <- Hs[((j-1)*d+1) : (j*d),]
    Hjs <- numeric(0)
    for (i in 1:nj) Hjs <- rbind(Hjs, Hj) # concatenate b/w matrices 
     
    fhat <- fhat + nj* dmvnorm.mixt(x=eval.points, mus=xj, Sigmas=Hjs,
                                    props=rep(1, nj)/nj)
  } 
  fhat <- fhat/n 
  
  return(list(x=x.pc, eval.points=eval.points, estimate=fhat, H=Hs))
}


###############################################################################
# Display kernel density estimate
#
# Parameters
# fhat - output from call to `kde'
###############################################################################

plot.kde <- function(x, display="slice", ...)
{ 
  fhat <- x
  d <- ncol(fhat$x)
  rm(x)

  if (d==2) 
    plotkde.2d(fhat, display=display, ...)
  else if (d== 3)
    plotkde.3d(fhat, display="rgl", ...)
  else 
    stop ("Plot function only available for 2 or 3-dimensional data")
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

plotkde.2d <- function(fhat, display="slice", cont=c(25,50,75), ncont=NULL,cex=0.7, 
    xlabs="x", ylabs="y", zlabs="Density function", theta=-30, phi=40, d=4,
    add=FALSE, drawlabels=TRUE, points.diff=TRUE, pch, ptcol="blue", lcol="black",
    ...)
{
  disp <- substr(display,1,1)
  if (!is.list(fhat$eval.points))
    stop("Need a grid of density estimates")

  if (is.data.frame(fhat$x))
  {
    xlabs <- names(fhat$x)[1]
    ylabs <- names(fhat$x)[2]
  }
  
  # perspective/wire-frame plot
  if (disp=="p")
    persp(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
          theta=theta, phi=phi, d=d, xlab=xlabs, ylab=ylabs, zlab=zlabs, ...)
  else
  {
    # slice/contour plot
    if (disp=="s")
    {
      # pre-clustered KDE
      if (is.list(fhat$x))
      {
        x.pc <- fhat$x
        d <- ncol(x.pc$x)
        nu <- length(x.pc$nclust)
        
        if (missing(pch)) pch <- 1:nu

        dobs <- kde(x.pc, fhat$H, eval.points=x.pc$x)$estimate
        hts <- quantile(dobs, prob = (100 - cont)/100)

        if (add)
          points(x.pc$x[,1], x.pc$x[,2], cex=cex, pch=pch[1], col=ptcol)
        else
          plot(x.pc$x[,1], x.pc$x[,2], type="n", xlab=xlabs, ylab=ylabs,
               col=ptcol,...)
      }
      # fixed bandwidth KDE
      else
      {
        if (missing(cex)) cex <- 0.7
        dobs <- kde(fhat$x, fhat$H, eval.points=fhat$x)$estimate 
        hts <- quantile(dobs, prob = (100 - cont)/100)

        if (add)
          points(fhat$x[,1], fhat$x[,2], cex=cex, col=ptcol)
        else
          plot(fhat$x[,1], fhat$x[,2], type="n", xlab=xlabs, ylab=ylabs, ...)
      }

      # compute and draw contours
      if (is.null(ncont))
        for (i in 1:length(cont)) 
        {
          scale <- cont[i]/hts[i]
          contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                  fhat$estimate*scale, level=hts[i]*scale, add=TRUE, 
                  drawlabels=drawlabels, col=lcol, ...)
        }
      else
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate,
                nlevel=ncont, add=TRUE, drawlabels=drawlabels, col=lcol, ...)
      
      # add points 
      if (is.list(fhat$x))
      {
        if (points.diff)
          for (j in 1:length(x.pc$nclust))
            if (is.vector(x.pc$x[x.pc$ind==j,]))
              points(x.pc$x[x.pc$ind==j,1], x.pc$x[x.pc$ind==j,2], cex=cex,
                     pch=pch[j], col=ptcol[j])
            else
              points(x.pc$x[x.pc$ind==j,], cex=cex, pch=pch[j], col=ptcol[j])
        else
          points(x.pc$x, cex=cex, pch=pch[1], col=ptcol[1])
      }
      else  
        points(fhat$x[,1], fhat$x[,2], cex=cex, col=ptcol)
    }
    # image plot
    else if (disp=="i")
      image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, 
            xlab=xlabs, ylab=ylabs, ...)
  }   
}



###############################################################################
# Display trivariate kernel density estimate
#
# Parameters 
# fhat - output from 'kde.grid'
# display - "persp" - perspective plot
#         - "slice" - contour plot
#         - "image" image plot
# cont - vector of contours to be plotted
###############################################################################


plotkde.3d <- function(fhat, display="rgl", cont=c(25,50,75), colors,
  alphavec, size=3, ptcol="blue", add=FALSE, origin=c(0,0,0),
  endpts, xlabs="x", ylabs="y", zlabs="z", drawpoints=TRUE, ...)

{
  dobs <- kde(fhat$x, fhat$H, eval.points=fhat$x)$estimate 
  hts <- quantile(dobs, prob = (100-cont)/100)
  nc <- length(cont)
  
  if (missing(colors))
    colors <- rev(heat.colors(nc))

  if (missing(endpts))
  {
    endpts <- rep(0,3)
    endpts[1] <-  max(fhat$eval.points[[1]])
    endpts[2] <-  max(fhat$eval.points[[2]])
    endpts[3] <-  max(fhat$eval.points[[3]])
  }

  #if (missing(xlim))
  #  xlim <- range(fhat$eval.points[[1]])
  #if (missing(ylim))
  #  ylim <- range(fhat$eval.points[[2]])
  #if (missing(zlim))
  #  zlim <- range(fhat$eval.points[[3]])

  #alph <- seq(alphalo, alphahi, length=nc)
  if (missing(alphavec))
    alphavec <- seq(0.1,0.5,length=nc)

  if (!add)
  {
    rgl.clear()
    rgl.bg(col="white")
    for (i in 1:nc) 
      contour3d(fhat$estimate, level=hts[nc-i+1], fhat$eval.points[[1]],
                fhat$eval.points[[2]], fhat$eval.points[[3]], add=(i>1),
                color=colors[i], alpha=alphavec[i], ...)
  }
  else
  {
    rgl.bg(col="white")
    for (i in 1:nc) 
      contour3d(fhat$estimate, level=hts[nc-i+1], fhat$eval.points[[1]],
                fhat$eval.points[[2]], fhat$eval.points[[3]], add=add,
                color=colors[i], alpha=alphavec[i], ...)
  }
   
  points3d.rh(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=ptcol)

  
  lines3d(c(origin[1],endpts[1]),rep(origin[2],2),rep(origin[3],2),size=3,
          color="black", add=TRUE)
  lines3d(rep(origin[1],2),c(origin[2],endpts[2]),rep(origin[3],2),size=3,
          color="black", add=TRUE)
  lines3d(rep(origin[1],2),rep(origin[2],2),c(origin[3],endpts[3]),size=3,
          color="black",add=TRUE)

  texts3d.rh(endpts[1]+0.1*abs(endpts[1]),origin[2],origin[3],xlabs,color="black",size=3)
  texts3d.rh(origin[1],endpts[2]+0.1*abs(endpts[2]),origin[3],ylabs,color="black",size=3)
  texts3d.rh(origin[1],origin[2],endpts[3]+0.1*abs(endpts[3]),zlabs,color="black",size=3)

  ### add labels for origin and axis limits
  #xlim.str <- toString(signif(xlim[2],3))
  #ylim.str <- toString(signif(ylim[2],3))
  #zlim.str <- toString(signif(zlim[2],3))
  #org.str <- paste("(",toString(signif(origin,3)),")", sep="")
   
  #texts3d.rh(xlim[2],origin[2],origin[3],xlim.str,color="black",size=3)
  #texts3d.rh(origin[1],ylim[2],origin[3],ylim.str,color="black",size=3)
  #texts3d.rh(origin[1],origin[2],zlim[2],zlim.str,color="black",size=3)
  #texts3d.rh(origin[1],origin[2],origin[2],org.str,color="black",size=3)
}


