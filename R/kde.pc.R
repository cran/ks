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
  gridx <- make.grid.ks(x, matrix.sqrt(Hmax), tol=supp, gridsize=gridsize) 
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
    xlabs, ylabs, zlabs="Density function", theta=-30, phi=40, d=4,
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

  x.names <- colnames(fhat$x) 
  if (!is.null(x.names))
  {
    if (missing(xlabs))
      xlabs <- x.names[1]
    if (missing(ylabs))
      ylabs <- x.names[2]
  }
  else
  {
    xlabs="x"
    ylabs="y"
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

