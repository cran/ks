###############################################################################
# Kernel discriminant analysis
###############################################################################


###############################################################################
# Find bandwidths for each class in training set
#
# Parameters
# x - data values
# group - group variable
# bw - type of bandwidth selector
# nstage, pilot, pre - parameters for plugin bandwidths
# diag - FALSE - use full b/w matrices
#      - TRUE - use diag b/w matrices
#
# Returns
# Matrix of bandwidths for each group in training set
###############################################################################

Hkda <- function(x, x.group, Hstart, bw="plugin", nstage=2, pilot="samse",
                 pre="sphere")
{
  d <- 2
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  bw <- substr(tolower(bw),1,1)
  Hs <- numeric(0)
  for (i in 1:m)
  {
    y <- x[x.group==grlab[i],]
    if (!missing(Hstart)) 
    {
      Hstarty <- Hstart[((i-1)*d+1) : (i*d),]
      if (bw=="l")        
          H <- Hlscv(y, Hstart=Hstarty)
      else if (bw=="s")
        H <- Hscv(y, pre=pre, Hstart=Hstarty)
      else if (bw=="p")
          H <- Hpi(y, nstage=nstage, pilot=pilot, pre=pre, Hstart=Hstarty)
        
    }
    else
    {
      if (bw=="l")
        H <- Hlscv(y)
      else if (bw=="s")
        H <- Hscv(y, pre=pre)
      else if (bw=="p")
        H <- Hpi(y, nstage=nstage, pilot=pilot, pre=pre)
    }
    Hs <- rbind(Hs, H)
  }

  return(Hs)   
}

Hkda.diag <- function(x, x.group, bw="plugin", nstage=2, pilot="samse",
                 pre="sphere")
{
  d <- 2
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  bw <- substr(tolower(bw),1,1)
  Hs <- numeric(0)
  for (i in 1:m)
  {
    y <- x[x.group==grlab[i],]
    if (bw=="l")
      H <- Hlscv.diag(y)
    else if (bw=="p") 
      H <- Hpi.diag(y, nstage=nstage, pilot=pilot, pre=pre)
    Hs <- rbind(Hs, H)
  }

  return(Hs)   
}



###############################################################################
# Classify data set according to discriminant analysis based on training data
#
# Parameter
# x - training data
# x.group - group variable for x
# y - data values to be classified
# Hs - bandwidth matrices
# prior.prob - prior probabilities
#
# Returns
# Group classification of data set y
###############################################################################

kda <- function(x, x.group, Hs, y, prior.prob=NULL)
{  
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.data.frame(y)) y <- as.matrix(y)
  gr <- sort(unique(x.group))

  # if prior.prob is NULL then use sample proportions
  if (is.null(prior.prob))
  {
    prior.prob <- rep(0, length(gr))
    for (j in 1:length(gr))
      prior.prob[j] <- length(which(x.group==gr[j]))
    prior.prob <- prior.prob/nrow(x)
  }
  
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

  ## Compute KDE and weighted KDE 
  m <- length(gr)
  fhat <- kda.kde(x, x.group, Hs, eval.points=y)
  fhat.wt <- matrix(0, ncol=m, nrow=nrow(y))  

  for (j in 1:m)
    fhat.wt[,j] <- fhat$est[[j]]* prior.prob[j]

  ## Assign y according largest weighted density value 
  disc.gr.temp <- apply(fhat.wt, 1, which.max)

  disc.gr <- gr
  for (j in 1:m)
  {
    ind <- which(disc.gr.temp==j)
    disc.gr[ind] <- gr[j]
  }
 
  return(disc.gr) 
}


###############################################################################
# Compares true group classification with an estimated one
#
# Parameters
# group - true group variable
# est.group - estimated group variable
#
# Returns
# List with components
# comp - cross-classification table of groupings - true groups are the rows,
#        estiamted groups are the columns
# error - total mis-classification rate
###############################################################################

compare <- function(x.group, est.group)
{
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  comp <- matrix(0, nr=m, nc=m)

  for (i in 1:m)
    for (j in 1:m)
      comp[i,j] <- sum((x.group==grlab[i]) & (est.group==grlab[j]))
  
  er <- 1 - sum(diag(comp))/sum(comp)
  colnames(comp) <- as.character(grlab)
  rownames(comp) <- as.character(grlab)
  comp <- cbind(comp, rowSums(comp))
  comp <- rbind(comp, colSums(comp))
  
  return(list(cross=comp, error=er))
 
}
                       

###############################################################################
# KDEs of individual densities for KDA
#
# Parameters
# x - data values
# group - group variable
# Hs - bandwidth matrices
#
# Returns
# List with components
# x - list of data values
# eval.points - evaluation points of dnesity estimate
# estimate - list of density estimate
# H - list of bandwidth matrices
##############################################################################

kda.kde <- function(x, x.group, Hs, gridsize=c(100,100), supp=3.7, eval.points=NULL)
{
  if (is.data.frame(x)) x <- as.matrix(x)
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  d <- ncol(x)
  
  # find largest bandwidth matrix to initialise grid
  detH <- vector() 
  for (j in 1:m)
    detH[j] <- det(Hs[((j-1)*d+1) : (j*d),])  
  Hmax.ind <- which.max(detH)
  Hmax <- Hs[((Hmax.ind-1)*d+1) : (Hmax.ind*d),]

  # initialise grid 
  gridx <- make.grid(x, matrix.sqrt(Hmax), tol=supp, gridsize=gridsize) 
  suppx <- make.supp(x, matrix.sqrt(Hmax), tol=supp)  
  grid.pts <- find.gridpts(gridx, suppx)
  fhat.list <- list()
  
  
  for (j in 1:m)
  {
    y <- x[x.group==grlab[j],]
    H <- Hs[((j-1)*d+1) : (j*d),]

    # compute individual density estimate
    if (is.null(eval.points))
    {
      y.grid.pts <- list()
      y.grid.pts$minx <- grid.pts$minx[x.group==grlab[j],]
      y.grid.pts$maxx <- grid.pts$maxx[x.group==grlab[j],]
      fhat.temp <- kde.grid.2d(y, H, gridx=gridx, supp=supp, grid.pts=y.grid.pts)
    }
    else
      fhat.temp <- kde.points(y, H, eval.points=eval.points)

    fhat.list$x <- c(fhat.list$x, list(y))
    fhat.list$eval.points <- fhat.temp$eval.points
    fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$est))
    fhat.list$H <- c(fhat.list$H, list(fhat.temp$H))
  }
  
  pr <- rep(0, length(grlab))
  for (j in 1:length(grlab))
    pr[j] <- length(which(x.group==grlab[j]))
  pr <- pr/nrow(x)
  fhat.list$prior.prob <- pr
  fhat.list$type <- "kernel"
  
  class(fhat.list) <- "dade"
  
  return(fhat.list)
}


##############################################################################
# Plot KDE of individual densities and partition
#
# Parameters
# fhat - output from `kda.kde'
# y - data points (separate from training data inside fhat)
# y.group - data group labels
# prior.prob - vector of prior porbabilities
# disp - "part" - plot partition
#      - "" - don't plot partition
#
##############################################################################


plot.dade <- function(x, y, y.group, prior.prob=NULL, display="part",
    cont=c(25,50,75), ncont=NULL, xlim, ylim, xlabs="x", ylabs="y",
    drawlabels=TRUE, cex=1, pch, lty, col, lcol, ...)
{ 
  fhat <- x
  rm(x)
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
  
  if (missing(y)) 
    plot(fhat$x[[1]], type="n", xlab=xlabs, ylab=ylabs, xlim=xlim, ylim=ylim, ...)
  else
    plot(y,  type="n", xlab=xlabs, ylab=ylabs, xlim=xlim, ylim=ylim, ...)


  if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob
  if (m != length(prior.prob))
    stop("Prior prob. vector not same length as number of components in fhat")
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

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
    
    if (missing(col)) col <- 2:(max(class.grid)+1)
    image(fhat$eval[[1]], fhat$eval[[2]], class.grid,col=col, xlim=xlim, ylim=ylim,
          add=TRUE, ...)
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
    if (missing(y))
      points(fhat$x[[j]], cex=cex, pch=pch[j])
    else 
    {
      if (missing(y.group))
        points(y, cex=cex)
      else
        points(y[y.group==y.gr[j],], cex=cex, pch=pch[j])
    
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

##############################################################################
# Plot density estimate of individual densities and classifcation region
# for parametric DA
#
# Parameters
# x - training data
# x.goup - training data group labels 
# y - data points 
# y.group - data group labels
# prior.prob - vector of prior porbabilities
# disp - "part" - plot partition
#      - "" - don't plot partition
##############################################################################


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


pda.pde <- function(x, x.group, gridsize=c(100,100), type="quad", xlim, ylim)
{
  n <- nrow(x)
  gr <- sort(unique(x.group))
  typ <- substr(type,1,1) 
  if (missing(xlim))
    xlim <- range(x[,1])
  if (missing(ylim))
    ylim <- range(x[,2])
  
  ex <- seq(xlim[1]-abs(xlim[1])/10, xlim[2]+abs(xlim[2])/10, length=gridsize[1])
  ey <- seq(ylim[1]-abs(ylim[1])/10, ylim[2]+abs(ylim[2])/10, length=gridsize[2])
  xy <- permute(list(ex, ey))

  fhat <- list()
  fhat$x <- list()
  fhat$eval.points <- list(ex, ey)
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
    
    if (typ=="q")
      dens <- dmvnorm.mixt(xy, mu=xbari, Sigma=Si, props=1)
    else if (typ=="l")
      dens <- dmvnorm.mixt(xy, mu=xbari, Sigma=S, props=1)
    fhat$estimate[[i]] <- matrix(dens, nc=length(ex), byrow=FALSE) 
  }

  prior.prob <- rep(0, length(gr))
  for (j in 1:length(gr))
    prior.prob[j] <- length(which(x.group==gr[j]))
  prior.prob <- prior.prob/nrow(x)
  fhat$prior.prob <- prior.prob
  
  fhat$type <- type
  
  class(fhat) <- "dade"
  
  return(fhat)
}

