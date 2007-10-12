###############################################################################
# Kernel discriminant analysis
###############################################################################


###############################################################################
# Find bandwidths for each class in training set, for 2- to 6-dim 
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
                 pre="sphere", binned=FALSE, bgridsize)
{
  d <- ncol(x)
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  bw <- substr(tolower(bw),1,1)
  Hs <- numeric(0)

  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)
  
  for (i in 1:m)
  {
    y <- x[x.group==grlab[i],]
    if (!missing(Hstart)) 
    {
      Hstarty <- Hstart[((i-1)*d+1) : (i*d),]
      if (bw=="l")        
        H <- Hlscv(y, Hstart=Hstarty)
      else if (bw=="s")
        H <- Hscv(y, pre=pre, Hstart=Hstarty, binned=binned, bgridsize=bgridsize)
      else if (bw=="p")
        H <- Hpi(y, nstage=nstage, pilot=pilot, pre=pre, Hstart=Hstarty,
                 binned=binned, bgridsize=bgridsize) 
    }
    else
    {
      if (bw=="l")
        H <- Hlscv(y)
      else if (bw=="s")
        H <- Hscv(y, pre=pre, binned=binned, bgridsize=bgridsize)
      else if (bw=="p")
        H <- Hpi(y, nstage=nstage, pilot=pilot, pre=pre, binned=binned, bgridsize=bgridsize)
    }
    Hs <- rbind(Hs, H)
  }

  return(Hs)   
}

Hkda.diag <- function(x, x.group, bw="plugin", nstage=2, pilot="samse",
                 pre="sphere", binned=FALSE, bgridsize)
{
  d <- ncol(x)
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  bw <- substr(tolower(bw),1,1)
  Hs <- numeric(0)

  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  for (i in 1:m)
  {
    y <- x[x.group==grlab[i],]
    if (bw=="l")
      H <- Hlscv.diag(y, binned=binned, bgridsize=bgridsize)
    else if (bw=="p") 
      H <- Hpi.diag(y, nstage=nstage, pilot=pilot, pre=pre, binned=binned, bgridsize=bgridsize)
    else if (bw=="s")
      H <- Hscv.diag(y, pre=pre, binned=binned, bgridsize=bgridsize)
    Hs <- rbind(Hs, H)
  }

  return(Hs)   
}



###############################################################################
# Classify data set according to discriminant analysis based on training data
# for 1- to 6-dim
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

kda <- function(x, x.group, Hs, hs, y, prior.prob=NULL)
{
  if (is.vector(x))
  {
    disc.gr <- kda.1d(x=x, x.group=x.group, hs=hs, y=y, prior.prob=prior.prob)
  }
  else
  {  
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.data.frame(y)) y <- as.matrix(y)
    gr <- sort(unique(x.group))

    d <- ncol(x)
   
    ## if prior.prob is NULL then use sample proportions
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
    fhat <- kda.kde(x, x.group, Hs=Hs, eval.points=y)
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
  }
  
  return(disc.gr) 
}

kda.1d <- function(x, x.group, hs, y, prior.prob=NULL)
{ 
  gr <- sort(unique(x.group))

  # if prior.prob is NULL then use sample proportions
  if (is.null(prior.prob))
  {
    prior.prob <- rep(0, length(gr))
    for (j in 1:length(gr))
      prior.prob[j] <- length(which(x.group==gr[j]))
    prior.prob <- prior.prob/length(x)
  }
  
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

  ## Compute KDE and weighted KDE 
  m <- length(gr)
  fhat <- kda.kde(x, x.group, hs=hs, eval.points=y)
  fhat.wt <- matrix(0, ncol=m, nrow=length(y))  
  
  for (j in 1:m)
    fhat.wt[,j] <- fhat$estimate[[j]]* prior.prob[j]

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
#        estimated groups are the columns
# error - total mis-classification rate
###############################################################################

compare <- function(x.group, est.group, by.group=FALSE)
{
  if (length(x.group)!=length(est.group))
    stop("Group label vectors not the same length")
  
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  comp <- matrix(0, nr=m, nc=m)
  
  for (i in 1:m)
    for (j in 1:m)
      comp[i,j] <- sum((x.group==grlab[i]) & (est.group==grlab[j]))  

  if (by.group)
  {
    er <- vector()
    for (i in 1:m)
      er[i] <- 1-comp[i,i]/rowSums(comp)[i]
    er <- matrix(er, nc=1)
    er <- rbind(er, 1 - sum(diag(comp))/sum(comp)) 
    rownames(er) <- c(as.character(paste(grlab, "(true)")), "Total")
    colnames(er) <- "error"
    
  }
  else 
    er <- 1 - sum(diag(comp))/sum(comp)
  
  comp <- cbind(comp, rowSums(comp))
  comp <- rbind(comp, colSums(comp))

  colnames(comp) <- c(as.character(paste(grlab, "(est.)")), "Total")
  rownames(comp) <- c(as.character(paste(grlab, "(true)")), "Total")

  return(list(cross=comp, error=er))
 
}

###############################################################################
# Computes cross-validated misclassification rates (for use when test data is
# not independent of training data) for KDA
#
# Parameters
# x - training data
# x.group - group variable for x
# y - data values to be classified
# Hs - bandwidth matrices
# prior.prob - prior probabilities
#
# Returns
# List with components
# comp - cross-classification table of groupings - true groups are the rows,
#        estimated groups are the columns
# error - total mis-classification rate
###############################################################################

compare.kda.cv <- function(x, x.group, bw="plugin",
    prior.prob=NULL, Hstart, by.group=FALSE, trace=FALSE, binned=FALSE,
    bgridsize, recompute=FALSE, ...)
{
  n <- nrow(x)
  d <- ncol(x)
  
  if (!missing(Hstart))
    H <- Hkda(x, x.group, bw=bw, Hstart=Hstart, binned=binned, bgridsize=bgridsize, ...)
  else
    H <- Hkda(x, x.group, bw=bw, binned=binned, bgridsize=bgridsize, ...)

  ### classify data x using KDA rules based on x itself
  kda.group <- kda(x, x.group, Hs=H, y=x, prior.prob=prior.prob)
  comp <- compare(x.group, kda.group)
 
  gr <- sort(unique(x.group)) 
  kda.cv.gr <- x.group
  
  for (i in 1:n)
  {
    H.mod <- H
    ### find group that x[i] belongs to 
    ind <- which(x.group[i]==gr)
    indx <- x.group==gr[ind]
    indx[i] <- FALSE

    if (recompute)
    {
      ## compute b/w matrix for that group with x[i] excluded
      if (!missing(Hstart))
      {  
        Hstart.temp <- Hstart[((ind-1)*d+1):(ind*d),]
        
        if (substr(bw,1,1)=="p")
          H.temp <- Hpi(x[indx,], Hstart=Hstart.temp, binned=binned, bgridsize=bgridsize,...)
        else if (substr(bw,1,1)=="s")
          H.temp <- Hscv(x[indx,],  Hstart=Hstart.temp, binned=binned, bgridsize=bgridsize,...)
        else if (substr(bw,1,1)=="l") 
          H.temp <- Hlscv(x[indx,],  Hstart=Hstart.temp, binned=binned, bgridsize=bgridsize,...)
      }
      else
      {
        if (substr(bw,1,1)=="p")
          H.temp <- Hpi(x[indx,], binned=binned, bgridsize=bgridsize, ...)
        else if (substr(bw,1,1)=="s")
          H.temp <- Hscv(x[indx,], binned=binned, bgridsize=bgridsize, ...)
        else if (substr(bw,1,1)=="l")
          H.temp <- Hlscv(x[indx,], binned=binned, bgridsize=bgridsize, ...) 
      }
      
      H.mod[((ind-1)*d+1):(ind*d),] <- H.temp
    }
    ## recompute KDA estimate of groups with x[i] excluded
      
    if (trace)
      cat(paste("Processing data item:", i, "\n"))

    kda.cv.gr[i] <- kda(x[-i,], x.group[-i], Hs=H.mod, y=x, prior.prob=prior.prob)[i]
  }
  
  return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}

###############################################################################
### Same as compare.kda.cv except uses diagonal b/w matrices
###############################################################################

compare.kda.diag.cv <- function(x, x.group, bw="plugin", prior.prob=NULL,
   by.group=FALSE, trace=FALSE, binned=FALSE, bgridsize, recompute=FALSE, ...)
{
  n <- nrow(x)
  d <- ncol(x)

  H <- Hkda.diag(x, x.group, bw=bw, binned=binned, bgridsize=bgridsize, ...)
  kda.group <- kda(x, x.group, Hs=H, y=x, prior.prob=prior.prob)
  comp <- compare(x.group, kda.group)
 
  gr <- sort(unique(x.group)) 
  kda.cv.gr <- x.group
  
  for (i in 1:n)
  {
    H.mod <- H

    if (recompute)
    {
      ind <- which(x.group[i]==gr)
      indx <- x.group==gr[ind]
      indx[i] <- FALSE
      if (substr(bw,1,1)=="p")
        H.temp <- Hpi.diag(x[indx,],  binned=binned, bgridsize=bgridsize, ...)
      else if (substr(bw,1,1)=="l")
        H.temp <- Hlscv.diag(x[indx,], binned=binned, bgridsize=bgridsize, ...)
      
      H.mod[((ind-1)*d+1):(ind*d),] <- H.temp
    }
    
    if (trace)
      cat(paste("Processing data item:", i, "\n"))

    kda.cv.gr[i] <- kda(x[-i,], x.group[-i], Hs=H.mod, y=x, prior.prob=prior.prob)[i]  
  }
  
  return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}



###############################################################################
# KDEs of individual densities for KDA - 1- to 3-dim
#
# Parameters
# x - data values
# group - group variable
# Hs - bandwidth matrices
#
# Returns
# List with components (class dade)
# x - list of data values
# eval.points - evaluation points of dnesity estimate
# estimate - list of density estimate
# H - list of bandwidth matrices
##############################################################################

kda.kde <- function(x, x.group, Hs, hs, prior.prob=NULL, gridsize, supp=3.7, eval.points=NULL, binned=FALSE, bgridsize)
{
  if (is.vector(x))
  {
    if (missing(gridsize))  gridsize <- 101
    fhat.list <- kda.kde.1d(x=x, x.group=x.group, hs=hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, eval.points=eval.points)
  }
  else
  {
    if (is.data.frame(x)) x <- as.matrix(x)
    grlab <- sort(unique(x.group))
    m <- length(grlab)
    d <- ncol(x)

    ## find largest bandwidth matrix to initialise grid
    detH <- vector() 
    for (j in 1:m)
      detH[j] <- det(Hs[((j-1)*d+1) : (j*d),])  
    Hmax.ind <- which.max(detH)
    Hmax <- Hs[((Hmax.ind-1)*d+1) : (Hmax.ind*d),]

    if (binned)
    {
      if (d > 4)
        stop("Binning only available for 1- to 4-dim data")
      
      if (missing(bgridsize))
        if (d==2)
          bgridsize <- rep(151,d)
        else if (d==3)
          bgridsize <- rep(51, d)
        else if (d==4)
          bgridsize <- rep(21, d)

      ## linear binning
      bin.par.x <- dfltCounts.ks(x, bgridsize, sqrt(diag(Hmax)), supp=supp)
      
      fhat.list <- list()
      
      for (j in 1:m)
      {
        xx <- x[x.group==grlab[j],]     
        H <- Hs[((j-1)*d+1) : (j*d),]     

        if (!identical(diag(diag(H)), H))
          stop("Binned estimation requires diagonal bandwidth matrices")
        
        bin.par <- dfltCounts.ks(xx, bgridsize, sqrt(diag(Hmax)), range.x=bin.par.x$range.x, supp=supp)
        
        if (is.null(eval.points))
        {
          fhat.temp <- drvkde.ks(x=bin.par$counts, drv=rep(0,d),
                                  bandwidth=sqrt(diag(H)),
                                  binned=TRUE, range.x=bin.par$range.x, se=FALSE)
          fhat.list$x <- c(fhat.list$x, list(xx))
          fhat.list$eval.points <- fhat.temp$x.grid
          fhat.temp$est[fhat.temp$est<0] <- 0
          fhat.temp$est <- zapsmall(fhat.temp$est, digits=20)
          fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$est))
          fhat.list$H <- c(fhat.list$H, list(H))
        }
        else
        {
          fhat.temp <- kde.points(xx, H, eval.points=eval.points)
          fhat.list$x <- c(fhat.list$x, list(xx))
          fhat.list$eval.points <- fhat.temp$eval.points
          fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$est))
          fhat.list$H <- c(fhat.list$H, list(fhat.temp$H))
        }
      }
    }
    else
    {
      if (missing(gridsize))
        if (d==2)
          gridsize <- rep(151,d)
        else if (d==3)
          gridsize <- rep(51, d)
        else if (d==4)
          gridsize <- rep(21, d)
        else
          gridsize <- rep(11, d)
       
      ## initialise grid 
      gridx <- make.grid.ks(x, matrix.sqrt(Hmax), tol=supp, gridsize=gridsize) 
      suppx <- make.supp(x, matrix.sqrt(Hmax), tol=supp)  
      grid.pts <- find.gridpts(gridx, suppx)
      fhat.list <- list()
      
      for (j in 1:m)
      {
        xx <- x[x.group==grlab[j],]
        H <- Hs[((j-1)*d+1) : (j*d),]
        
        ## compute individual density estimate
        if (is.null(eval.points))
        {
          xx.grid.pts <- list()
          xx.grid.pts$minx <- grid.pts$minx[x.group==grlab[j],]
          xx.grid.pts$maxx <- grid.pts$maxx[x.group==grlab[j],]
          if (d==2)
            fhat.temp <- kde.grid.2d(xx, H, gridx=gridx, supp=supp, grid.pts=xx.grid.pts)
          else if (d==3)
            fhat.temp <- kde.grid.3d(xx, H, gridx=gridx, supp=supp, grid.pts=xx.grid.pts)
        }
        else
          fhat.temp <- kde.points(xx, H, eval.points=eval.points)
        
        fhat.list$x <- c(fhat.list$x, list(xx))
        fhat.list$eval.points <- fhat.temp$eval.points
        fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$est))
        fhat.list$H <- c(fhat.list$H, list(fhat.temp$H))
      }
    }
    
    fhat.list$x.group <- x.group
    pr <- rep(0, length(grlab))
    for (j in 1:length(grlab))
      pr[j] <- length(which(x.group==grlab[j]))
    pr <- pr/nrow(x)
    fhat.list$prior.prob <- pr

    fhat.list$binned <- binned
    
    class(fhat.list) <- "kda.kde"
  }
  
  return(fhat.list)
}

kda.kde.1d <- function(x, x.group, hs, prior.prob, gridsize, supp, eval.points)
{
  grlab <- sort(unique(x.group))
  m <- length(grlab)

  hmax <- max(hs)
  range.x <- c(min(x) - supp*hmax, max(x) + supp*hmax)
  fhat.list <- list()
  for (j in 1:m)
  {
    xx <- x[x.group==grlab[j]]
    h <- hs[j]

    ## compute individual density estimate
    if (is.null(eval.points))
    {
      fhat.temp <- bkde(x=xx, bandwidth=h, gridsize=gridsize, range.x=range.x)
      fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$y))
    }
    else
    {
      fhat.temp <- kde.points.1d(x=xx, h=h, eval.points=eval.points)
      fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$estimate))
    }
    fhat.list$x <- c(fhat.list$x, list(xx))
    fhat.list$eval.points <- fhat.temp$x
    
    fhat.list$h <- c(fhat.list$h, h)
  }

  fhat.list$x.group <- x.group
  if (is.null(prior.prob))
  {
    pr <- rep(0, length(grlab))
    for (j in 1:length(grlab))
      pr[j] <- length(which(x.group==grlab[j]))
    pr <- pr/length(x)
    fhat.list$prior.prob <- pr
  }
  else
    fhat.list$prior.prob <- prior.prob
  #fhat.list$type <- "kernel"
  
  class(fhat.list) <- "kda.kde"
  
  return(fhat.list)
  
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


plot.kda.kde <- function(x, y, y.group, ...) 
{
  if (is.vector(x$x[[1]]))
    plotkda.kde.1d(x=x, y=y, y.group=y.group, ...)
  else
  {  
    d <- ncol(x$x[[1]])
    
    if (d==2)
      plotkda.kde.2d(x=x, y=y, y.group=y.group, ...) 
    else if (d==3)
      ##warning("RGL 3-d plotting temporarily disabled")  
       plotkda.kde.3d(x=x, y=y, y.group=y.group, ...) 
  }
}


plotkda.kde.1d <- function(x, y, y.group, prior.prob=NULL, xlim, ylim, xlab="x", ylab="Weighted density function", drawpoints=TRUE, lty, lcol, col, ptcol, ...)
{ 
  fhat <- x
  
  m <- length(fhat$x)
  ##type <- substr(fhat$type,1,1)
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
  {
    rug(fhat$eval.points[eval.points.gr==j], col=col[j])
    if (drawpoints)
      rug(fhat$x[[j]], col=ptcol[j], ticksize=-0.03)
  }
  
  if (!missing(y.group)) y.gr <- sort(unique(y.group))
  if (!missing(y))
    for (j in 1:length(y.gr))
      rug(y[y.group==y.gr[j]], col=ptcol[j], ticksize=-0.03)
}


plotkda.kde.2d <- function(x, y, y.group, prior.prob=NULL, 
    cont=c(25,50,75), abs.cont, xlim, ylim, xlab, ylab,
    drawpoints=FALSE, drawlabels=TRUE,  cex=1, pch, lty, col, lcol, ptcol, ...)
{ 
  fhat <- x
  
  d <- 2
  m <- length(fhat$x)
  ##type <- substr(fhat$type,1,1)
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
  box()
  
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

  if (fhat$binned)
    bin.par.xx <- dfltCounts.ks(xx, gridsize=dim(fhat$est[[j]]), sqrt(diag(fhat$H[[j]])), supp=3.7)
  
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
    if (missing(abs.cont))
    {
      if (fhat$binned)
      {
        bin.par <- dfltCounts.ks(fhat$x[[j]], gridsize=dim(fhat$est[[j]]), sqrt(diag(fhat$H[[j]])), supp=3.7)
        
        dobs <- c(dobs, rep(fhat$estimate[[j]], round(bin.par.xx$counts,0))* prior.prob[j])
      }
      else
        dobs <- c(dobs, kde(x=fhat$x[[j]], H=fhat$H[[j]], eval.points=xx)$estimate*prior.prob[j])
    }
  }

  if (missing(abs.cont))
    hts <- quantile(dobs, prob = (100 - cont)/100)
  else
    hts <- abs.cont
  
  for (i in 1:length(hts)) 
  {
    scale <- cont[i]/hts[i]
    if (missing(abs.cont))
      for (j in 1:m)
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                fhat$estimate[[j]]*scale, level=hts[i]*scale, add=TRUE, 
                drawlabels=drawlabels, lty=lty[j], col=lcol[j], ...)
    else
      for (j in 1:m)
        contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                fhat$estimate[[j]], level=hts[i], add=TRUE, 
                drawlabels=drawlabels, lty=lty[j], col=lcol[j], ...)
  }
}



plotkda.kde.3d <- function(x, y, y.group, prior.prob=NULL,
    cont=c(25,50), abs.cont, colors, alphavec, xlab, ylab, zlab,
    drawpoints=FALSE, size=3, ptcol="blue", ...)
{   
  fhat <- x
   
  d <- 3
  m <- length(fhat$x)
 
  eval1 <- fhat$eval.points[[1]]
  eval2 <- fhat$eval.points[[2]]
  eval3 <- fhat$eval.points[[3]]

  if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob
  if (m != length(prior.prob))
    stop("Prior prob. vector not same length as number of components in fhat")
  if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
    stop("Sum of prior weights not equal to 1")

  ##if (missing(endpts))
  ##{
  ##  endpts <- rep(0,3)
  ##  endpts[1] <-  max(fhat$eval.points[[1]])
  ##  endpts[2] <-  max(fhat$eval.points[[2]])
  ##  endpts[3] <-  max(fhat$eval.points[[3]])
  ##}

  x.names <- colnames(fhat$x[[1]])

  if (missing(xlab))
    if (is.null(x.names)) xlab <- "x" else xlab <- x.names[1]
  if (missing(ylab))
    if (is.null(x.names)) ylab <- "y" else ylab <- x.names[2]
  if (missing(zlab))
    if (is.null(x.names)) zlab <- "z" else zlab <- x.names[3]
             
  dobs <- numeric(0)
  xx <- numeric(0)

  for (j in 1:m)
    xx <- rbind(xx, fhat$x[[j]])

  x.gr <- sort(unique(fhat$x.group))

  if (fhat$binned)
    bin.par.xx <- dfltCounts.ks(xx, gridsize=dim(fhat$est[[j]]), sqrt(diag(fhat$H[[j]])), supp=3.7)

  if (missing(abs.cont))
  {
    for (j in 1:m)
      if (fhat$binned)
      {
        bin.par <- dfltCounts.ks(fhat$x[[j]], gridsize=dim(fhat$est[[j]]), sqrt(diag(fhat$H[[j]])), supp=3.7)
        
        dobs <- c(dobs, rep(fhat$estimate[[j]], round(bin.par.xx$counts,0))* prior.prob[j])
      }
      else
        dobs <- c(dobs, kde(x=fhat$x[[j]], H=fhat$H[[j]], eval.points=xx)$estimate*prior.prob[j])
    
    hts <- quantile(dobs, prob = (100-cont)/100)
  }
  else
    hts <- abs.cont

  ncont <- length(hts)
  
  if (missing(alphavec)) alphavec <- seq(0.1,0.3,length=ncont)
  if (missing(colors)) colors <- rainbow(m)
  if (missing(ptcol)) ptcol <- rep("blue", m)
  if (length(ptcol)==1) ptcol <- rep(ptcol, m)
  
  clear3d()
  bg3d(color="white")

  plot3d(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]],
         z=fhat$eval.points[[3]], type="n", xlab=xlab, ylab=ylab, zlab=zlab,
         ...)
  
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
  
  ##lines3d(c(origin[1],endpts[1]),rep(origin[2],2),rep(origin[3],2),size=3,
  ##        color="black", alpha=1)
  ##lines3d(rep(origin[1],2),c(origin[2],endpts[2]),rep(origin[3],2),size=3,
  ##        color="black", alpha=1)
  ##lines3d(rep(origin[1],2),rep(origin[2],2),c(origin[3],endpts[3]),size=3,
  ##        color="black", alpha=1)

  ##texts3d(endpts[1],origin[2],origin[3],xlab,color="black",size=3, alpha=1)
  ##texts3d(origin[1],endpts[2],origin[3],ylab,color="black",size=3, alpha=1)
  ##texts3d(origin[1],origin[2],endpts[3],zlab,color="black",size=3, alpha=1)
}



