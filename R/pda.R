###############################################################################
# Computes cross-validated misclassification rates (for use when test data =
# training data) for parametric DA
#
# Parameter
# x - training data
# x.group - group variable for x
# prior.prob - prior probabilities
# type - "line" - linear disc.
#      - "quad" - quadratic disc.
###############################################################################

compare.pda.cv <- function(x, x.group, type="quad", prior.prob=NULL,
    by.group=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)

  pda.group <- pda(x, x.group, x, prior.prob=prior.prob, type=type)
  comp <- compare(x.group, pda.group)
  
  pda.cv.gr <- vector()
  for (i in 1:n)
    pda.cv.gr[i] <-
      as.vector(pda(x[-i,], x.group[-i], x, prior.prob=prior.prob, type=type))[i]

  return(compare(x.group, pda.cv.gr, by.group=by.group)) 
}
##############################################################################
# Plot KDE of individual densities and partition - only for 2-dim
#
# Parameters
# fhat - output from `kda.kde'
# y - data points (separate from training data inside fhat)
# y.group - data group labels
# prior.prob - vector of prior probabilities
# disp - "part" - plot partition
#      - "" - don't plot partition
##############################################################################


plot.dade <- function(x, y, y.group, ...) #prior.prob=NULL, display="part",
    #cont=NULL, ncont=NULL, ...)
{
  if (is.vector(x$x[[1]]))
    plotdade.1d(x=x, y=y, y.group=y.group, ...)
  else
  {  
    d <- ncol(x$x[[1]])
    
    if (d==2)
      plotdade.2d(x=x, y=y, y.group=y.group, ...) 
    else if (d==3)
      plotdade.3d(x=x, y=y, y.group=y.group, ...) 
  }
}


plotdade.1d <- function(x, y, y.group, prior.prob=NULL, xlim, ylim, xlab="x", ylab="Weighted density function", drawpoints=TRUE, lty, lcol, col, ptcol, ...)
{ 
  fhat <- x
  
  m <- length(fhat$x)
  type <- substr(fhat$type,1,1)
  eval1 <- fhat$eval.points

  if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob
  
  if (m != length(prior.prob))
    stop("Prior prob. vector not same length as number of components in fhat")
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

  weighted.fhat <- matrix(0, nrow=length(fhat$eval.points), ncol=m) 
  for (j in 1:m)
    weighted.fhat[,j] <- fhat$estimate[[j]]*fhat$prior.prob[j]
  
  if (missing(xlim)) xlim <- range(fhat$eval.points)
  if (missing(ylim)) ylim <- range(weighted.fhat)
  if (missing(lty)) lty <- 1:m
  if (missing(lcol)) lcol <- 1:m
  if (missing(col)) col <- 1:m
  if (missing(ptcol)) ptcol <- rep("blue", m)
  

  plot(fhat$eval.points, weighted.fhat[,1], type="l", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, lty=lty[1], col=lcol[1], ...)
  
  if (m > 1)
    for (j in 2:m)
      lines(fhat$eval.points, weighted.fhat[,j], lty=lty[j], col=lcol[j], ...)

  eval.points.gr <- apply(weighted.fhat, 1, which.max)
  for (j in 1:m)
    rug(fhat$eval.points[eval.points.gr==j], col=col[j])
 
  if (!missing(y.group)) y.gr <- sort(unique(y.group))
  if (!missing(y))
    for (j in 1:length(y.gr))
      rug(y[y.group==y.gr[j]], col=ptcol[j], ticksize=-0.03)
}


plotdade.2d <- function(x, y, y.group, prior.prob=NULL, display="part",
    cont=c(25,50,75), ncont=NULL, xlim, ylim, xlab, ylab,
    drawpoints=TRUE, drawlabels=TRUE, cex=1, pch, lty, col, lcol, ptcol, ...)
{ 
  fhat <- x
  
  d <- 2
  m <- length(fhat$x)
  type <- substr(fhat$type,1,1)
  eval1 <- fhat$eval.points[[1]]
  eval2 <- fhat$eval.points[[2]]
  
  if (missing(xlim)) xlim <- c(min(eval1), max(eval1))
  if (missing(ylim)) ylim <- c(min(eval2), max(eval2))
  if (missing(pch)) pch <- 1:m
  if (missing(lty)) lty <- 1:m
  if (missing(lcol)) lcol <- rep(1, m)
  if (missing(ptcol)) ptcol <- rep("blue", m)

  x.names <- colnames(fhat$x[[1]]) 
  if (!is.null(x.names))
  {
    if (missing(xlab)) xlab <- x.names[1]
    if (missing(ylab)) ylab <- x.names[2]
  }
  else
  {
    xlab="x"
    ylab="y"
  }

  if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob

  if (m != length(prior.prob))
    stop("Prior prob. vector not same length as number of components in fhat")
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")
  
  if (missing(y)) 
    plot(fhat$x[[1]], type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  else
    plot(y, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  
  if (display=="part")
  {
    class.grid <- array(0, dim=dim(fhat$est[[1]]))
    temp <- matrix(0, ncol=length(fhat$est), nrow=nrow(fhat$est[[1]]))
    for (j in 1:ncol(fhat$est[[1]]))
    {
      for (k in 1:length(fhat$est))
        temp[,k] <- fhat$est[[k]][,j]* prior.prob[k]
      class.grid[,j] <- max.col(temp)

    }
    
    if (missing(col)) col <- heat.colors(m)
    image(fhat$eval[[1]], fhat$eval[[2]], class.grid,col=col, xlim=xlim,
          ylim=ylim, add=TRUE, ...)
  }
  
  dobs <- numeric(0)
  xx <- numeric(0)
  S <- 0
  n <- 0 
  for (j in 1:m)
  {
    n <- nrow(fhat$x[[j]])
    S <- S + nrow(fhat$x[[j]]) * var(fhat$x[[j]])
  }
  S <- S/(n - m)


  for (j in 1:m)
    xx <- rbind(xx, fhat$x[[j]])
  if (!missing(y.group)) y.gr <- sort(unique(y.group))

  for (j in 1:m)
  {
    if (drawpoints)
    {
      if (missing(y))
        points(fhat$x[[j]], cex=cex, pch=pch[j], col=ptcol[1])
      else 
      {
        if (missing(y.group))
          points(y, cex=cex, col=ptcol[1])
        else
          points(y[y.group==y.gr[j],], cex=cex, pch=pch[j], col=ptcol[j]) 
      }
    }
    
    if (type=="k")
      dobs <- c(dobs, kde(fhat$x[[j]], fhat$H[[j]], eval.points=xx)$estimate
                * prior.prob[j])
    else if (type=="q")
    {
      xbarj <- apply(fhat$x[[j]], 2, mean)
      Sj <- var(fhat$x[[j]])
      dobs <- c(dobs, dmvnorm.mixt(x=fhat$x[[j]], mus=xbarj, Sigmas=Sj, props=1)
                * prior.prob[j])
    }
    else if (type=="l")
    {
      xbarj <- apply(fhat$x[[j]], 2, mean)      
      dobs <- c(dobs, dmvnorm.mixt(x=fhat$x[[j]], mus=xbarj, Sigmas=S, props=1)
                * prior.prob[j])
    }
  }

  hts <- quantile(dobs, prob = (100 - cont)/100)
  
  if (is.null(ncont))
    for (i in 1:length(cont)) 
    {
      scale <- cont[i]/hts[i]
      for (j in 1:m)
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                fhat$estimate[[j]]*scale, level=hts[i]*scale, add=TRUE, 
                drawlabels=drawlabels, lty=lty[j], col=lcol[j], ...)
    }
  else
    for (j in 1:m)
      contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
              fhat$estimate[[j]], add=TRUE,  drawlabels=drawlabels,
              nlevel=ncont, lty=lty[j], col=lcol[j], ...)
  
}



plotdade.3d <- function(x, y, y.group, prior.prob=NULL, display="rgl",
    cont=c(25,50), colors, alphavec, origin=c(0,0,0),
    endpts, xlab, ylab, zlab, drawpoints=TRUE, size=3,
    ptcol, ...)
{ 
  fhat <- x
   
  d <- 3
  m <- length(fhat$x)
  type <- substr(fhat$type,1,1)
  eval1 <- fhat$eval.points[[1]]
  eval2 <- fhat$eval.points[[2]]
  eval3 <- fhat$eval.points[[3]]

  if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob
  if (m != length(prior.prob))
    stop("Prior prob. vector not same length as number of components in fhat")
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

  if (missing(endpts))
  {
    endpts <- rep(0,3)
    endpts[1] <-  max(fhat$eval.points[[1]])
    endpts[2] <-  max(fhat$eval.points[[2]])
    endpts[3] <-  max(fhat$eval.points[[3]])
  }
  x.names <- colnames(fhat$x[[1]]) 
  if (!is.null(x.names))
  {
    if (missing(xlab)) xlab <- x.names[1]
    if (missing(ylab)) ylab <- x.names[2]
    if (missing(zlab)) zlab <- x.names[3]
  }
  else
  {
    xlab="x"
    ylab="y"
    zlab="z"
  }
  
  ncont <- length(cont)
  
  if (missing(alphavec)) alphavec <- seq(0.1,0.3,length=ncont)
  if (missing(colors)) colors <- heat.colors(m)
  if (missing(ptcol)) ptcol <- rep("blue", m)
                 
  dobs <- numeric(0)
  xx <- numeric(0)
  S <- 0
  n <- 0 

  for (j in 1:m)
  {
    n <- nrow(fhat$x[[j]])
    S <- S + nrow(fhat$x[[j]]) * var(fhat$x[[j]])
  }
  S <- S/(n - m)


  for (j in 1:m)
    xx <- rbind(xx, fhat$x[[j]])

  #if (!missing(y.group))
  x.gr <- sort(unique(fhat$x.group))

  for (j in 1:m)
  {
    if (type=="k")
      dobs <- c(dobs, kde(fhat$x[[j]], fhat$H[[j]], eval.points=xx)$estimate
                * prior.prob[j])
    else if (type=="q")
    {
      xbarj <- apply(fhat$x[[j]], 2, mean)
      Sj <- var(fhat$x[[j]])
      dobs <- c(dobs, dmvnorm.mixt(x=fhat$x[[j]], mus=xbarj, Sigmas=Sj, props=1)
                * prior.prob[j])
    }
    else if (type=="l")
    {
      xbarj <- apply(fhat$x[[j]], 2, mean)      
      dobs <- c(dobs, dmvnorm.mixt(x=fhat$x[[j]], mus=xbarj, Sigmas=S, props=1)
                * prior.prob[j])
    }
  }
  
  hts <- quantile(dobs, prob = (100-cont)/100)

  clear3d()
  bg3d(color="white")
  
  for (j in 1:m)
  {
    for (i in 1:ncont) 
      contour3d(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]],
                z=fhat$eval.points[[3]], f=fhat$estimate[[j]],
                level=hts[ncont-i+1],
                add=TRUE, alpha=alphavec[i], color=colors[j],...)

    if (drawpoints)   ## plot points
    {
      if (missing(y))
        points3d(fhat$x[[j]][,1], fhat$x[[j]][,2], fhat$x[[j]][,3],
                    color=ptcol[j], size=size, alpha=1)
      else
      {
        if (missing(y.group))
          points3d(y[,1], y[,2], y[,3], color=ptcol, size=size, alpha=1)
        else
        {
          y.temp <- y[y.group==x.gr[j],]
          if (nrow(y.temp)>0)
            points3d(y.temp[,1], y.temp[,2], y.temp[,3], color=ptcol[j], size=size, alpha=1)
        }
      }
    }
  }
  
  lines3d(c(origin[1],endpts[1]),rep(origin[2],2),rep(origin[3],2),size=3,
          color="black", alpha=1)
  lines3d(rep(origin[1],2),c(origin[2],endpts[2]),rep(origin[3],2),size=3,
          color="black", alpha=1)
  lines3d(rep(origin[1],2),rep(origin[2],2),c(origin[3],endpts[3]),size=3,
          color="black", alpha=1)

  texts3d(endpts[1],origin[2],origin[3],xlab,color="black",size=3, alpha=1)
  texts3d(origin[1],endpts[2],origin[3],ylab,color="black",size=3, alpha=1)
  texts3d(origin[1],origin[2],endpts[3],zlab,color="black",size=3, alpha=1)
}



###############################################################################
# Classify data set according to discriminant analysis based on training data
# using linear or quadratic discrimination 
#
# Parameter
# x - training data
# x.group - group variable for x
# y - data values to be classified
# prior.prob - prior probabilities
# type - "line" - linear disc.
#      - "quad" - quadratic disc.
#
# Returns
# Group classification of data set y
###############################################################################


pda <- function(x, x.group, y, prior.prob=NULL, type="quad")
{  
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(y)) y <- as.matrix(y)
  gr <- sort(unique(x.group))
  
  if (is.null(prior.prob))
  {
    prior.prob <- rep(0, length(gr))
    for (j in 1:length(gr))
      prior.prob[j] <- length(which(x.group==gr[j]))
    prior.prob <- prior.prob/nrow(x)
  }
  
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")
  
  if (substr(type,1,1)=="q")
  {
    x.qda <- qda(x, x.group)
    disc.gr <- predict(x.qda, y, dimen=ncol(x), prior=prior.prob)$class
  }
  else
  {
    x.lda <- lda(x, x.group)
    disc.gr <- predict(x.lda, y, dimen=ncol(x), prior=prior.prob)$class 
  }
  return(disc.gr) 
}

##############################################################################
# Compute density estimate of individual densities and classifcation region
# for parametric DA - for 2-d and 3-d only
#
# Parameters
# x - training data
# x.group - training data group labels 
# type - "line" - linear disc.
#      - "quad" - quadratic disc.
#
# Returns 
# List with components (class dade)
# x - list of data values
# eval.points - evaluation points of dnesity estimate
# estimate - list of density estimate 
##############################################################################

pda.pde <- function(x, x.group, gridsize, type="quad", xlim, ylim, zlim)
{
  n <- nrow(x)
  d <- ncol(x)
  
  gr <- sort(unique(x.group))
  typ <- substr(type,1,1) 
  if (missing(xlim)) xlim <- range(x[,1])
  if (missing(ylim)) ylim <- range(x[,2])
  if (missing(zlim) & d==3) zlim <- range(x[,3])
  if (missing(gridsize)) gridsize <- rep(100,d)

  if (d==2)
  {  
    ex <- seq(xlim[1]-abs(xlim[1])/10, xlim[2]+abs(xlim[2])/10, length=gridsize[1])
    ey <- seq(ylim[1]-abs(ylim[1])/10, ylim[2]+abs(ylim[2])/10, length=gridsize[2])
    xy <- permute(list(ex, ey))
  }
  else if (d==3)
  {
    ex <- seq(xlim[1]-abs(xlim[1])/10, xlim[2]+abs(xlim[2])/10, length=gridsize[1])
    ey <- seq(ylim[1]-abs(ylim[1])/10, ylim[2]+abs(ylim[2])/10, length=gridsize[2])
    ez <- seq(zlim[1]-abs(zlim[1])/10, zlim[2]+abs(ylim[2])/10, length=gridsize[3])
    xyz <- permute(list(ex, ey, ez))
  }
  
  fhat <- list()
  fhat$x <- list()
  if (d==2)
    fhat$eval.points <- list(ex, ey)
  else if(d==3)
    fhat$eval.points <- list(ex, ey, ez)

  fhat$estimate <- list()
  dens.mat <- list()

 
  S <- 0
  for (i in 1:length(gr))
  {
    xi <- x[x.group==gr[i],]
    Si <- var(xi)
    S <- S + nrow(xi)*Si
  }
  S <- S/(n - length(gr))
 
  for (i in 1:length(gr))
  {
    xi <- x[x.group==gr[i],]
    xbari <- apply(xi, 2, mean)
    Si <- var(xi)
    fhat$x[[i]] <- xi

    if (d==2)
    {
      if (typ=="q")
        dens <- dmvnorm.mixt(xy, mu=xbari, Sigma=Si, props=1)
      else if (typ=="l")
        dens <- dmvnorm.mixt(xy, mu=xbari, Sigma=S, props=1)
      fhat$estimate[[i]] <- matrix(dens, nc=gridsize[1], byrow=FALSE)
    }
    else if (d==3)
    {
      if (typ=="q")
        dens <- dmvnorm.mixt(xyz, mu=xbari, Sigma=Si, props=1)
      else if (typ=="l")
        dens <- dmvnorm.mixt(xyz, mu=xbari, Sigma=S, props=1)
      fhat$estimate[[i]] <- array(dens, dim=gridsize)
    }
  }

  fhat$x.group <- x.group
  prior.prob <- rep(0, length(gr))
  for (j in 1:length(gr))
    prior.prob[j] <- length(which(x.group==gr[j]))
  prior.prob <- prior.prob/nrow(x)
  fhat$prior.prob <- prior.prob
  
  fhat$type <- type
  
  class(fhat) <- "dade"
  
  return(fhat)
}




