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
                 pre="sphere")
{
  d <- ncol(x)
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
  d <- ncol(x)
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
# for 2- to 6-dim
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
#        estimated groups are the columns
# error - total mis-classification rate
###############################################################################

compare <- function(x.group, est.group, by.group=FALSE)
{
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
    prior.prob=NULL, Hstart, by.group=FALSE, trace=FALSE, ...)
{
  n <- nrow(x)
  d <- ncol(x)
  
  if (!missing(Hstart))
    H <- Hkda(x, x.group, bw=bw, Hstart=Hstart, ...)
  else
    H <- Hkda(x, x.group, bw=bw, ...)

  ### classify data x using KDA rules based on x itself
  kda.group <- kda(x, x.group, H, x, prior.prob=prior.prob)
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

    ### compute b/w matrix for that group with x[i] excluded
    if (!missing(Hstart))
    {  
      Hstart.temp <- Hstart[((ind-1)*d+1):(ind*d),]
      
      if (substr(bw,1,1)=="p")
        H.temp <- Hpi(x[indx,], Hstart=Hstart.temp, ...)
      else if (substr(bw,1,1)=="s")
        H.temp <- Hscv(x[indx,],  Hstart=Hstart.temp,...)
      else if (substr(bw,1,1)=="l")
        H.temp <- Hlscv(x[indx,],  Hstart=Hstart.temp,...)
    }
    else
    {
      if (substr(bw,1,1)=="p")
        H.temp <- Hpi(x[indx,],  ...)
      else if (substr(bw,1,1)=="s")
        H.temp <- Hscv(x[indx,], ...)
      else if (substr(bw,1,1)=="l")
        H.temp <- Hlscv(x[indx,], ...)
    }
      
    H.mod[((ind-1)*d+1):(ind*d),] <- H.temp

    ### recompute KDA estimate of groups with x[i] excluded

    if (trace)
      cat(paste("Processing data item:", i, "\n"))
    kda.cv.gr[i] <- kda(x[-i,], x.group[-i], H.mod, x, prior.prob=prior.prob)[i]
  }
  
  return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}

###############################################################################
### Same as compare.kda.cv except uses diagonal b/w matrices
###############################################################################

compare.kda.diag.cv <- function(x, x.group, bw="plugin", prior.prob=NULL,
   by.group=FALSE, trace=FALSE,...)
{
  n <- nrow(x)
  d <- ncol(x)

  H <- Hkda.diag(x, x.group, bw=bw, ...)
  kda.group <- kda(x, x.group, H, x, prior.prob=prior.prob)
  comp <- compare(x.group, kda.group)
 
  gr <- sort(unique(x.group)) 
  kda.cv.gr <- x.group
  
  for (i in 1:n)
  {
    H.mod <- H

    ind <- which(x.group[i]==gr)
    indx <- x.group==gr[ind]
    indx[i] <- FALSE
    if (substr(bw,1,1)=="p")
      H.temp <- Hpi.diag(x[indx,],  ...)
    else if (substr(bw,1,1)=="l")
      H.temp <- Hlscv.diag(x[indx,], ...)
    
    #H[((ind-1)*d+1):(ind*d),] <- H.temp
    H.mod[((ind-1)*d+1):(ind*d),] <- H.temp
    if (trace)
      cat(paste("Processing data item:", i, "\n"))
    kda.cv.gr[i] <- kda(x[-i,], x.group[-i], H.mod, x, prior.prob=prior.prob)[i]  
  }
  return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}

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


###############################################################################
# KDEs of individual densities for KDA - only for 2-dim 
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

kda.kde <- function(x, x.group, Hs, gridsize, supp=3.7, eval.points=NULL)
{
  if (is.data.frame(x)) x <- as.matrix(x)
  grlab <- sort(unique(x.group))
  m <- length(grlab)
  d <- ncol(x)
   
  if (missing(gridsize)) 
    gridsize <- rep(100,d)
   
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
    xx <- x[x.group==grlab[j],]
    H <- Hs[((j-1)*d+1) : (j*d),]

    # compute individual density estimate
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
# Plot KDE of individual densities and partition - only for 2-dim
#
# Parameters
# fhat - output from `kda.kde'
# y - data points (separate from training data inside fhat)
# y.group - data group labels
# prior.prob - vector of prior porbabilities
# disp - "part" - plot partition
#      - "" - don't plot partition
##############################################################################


plot.dade <- function(x, y, y.group, prior.prob=NULL, display="part",
    cont=c(25,50,75), ncont=NULL, ...)
{
  d <- ncol(x$x[[1]])

  if (d==2)
    plotdade.2d(x, y, y.group, prior.prob=prior.prob, display=display,
                cont=cont, ncont=ncont, ...)
  else if (d==3)
    plotdade.3d(x, y, y.group, prior.prob=prior.prob, display="rgl",
                cont=cont, ...)
}

plotdade.2d <- function(x, y, y.group, prior.prob=NULL, display="part",
    cont=c(25,50,75), ncont=NULL, xlim, ylim, xlabs="x", ylabs="y",
    drawlabels=TRUE, cex=1, pch, lty, col, lcol, ptcol="blue", ...)
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
    if (missing(y))
      points(fhat$x[[j]], cex=cex, pch=pch[j], col=ptcol[1])
    else 
    {
      if (missing(y.group))
        points(y, cex=cex, col=ptcol[1])
      else
        points(y[y.group==y.gr[j],], cex=cex, pch=pch[j], col=ptcol[j])
    
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
    cont=c(25,50),  colors, alphavec, origin=c(0,0,0),
    endpts, xlabs="x", ylabs="y", zlabs="z", drawpoints=TRUE, size=3,
    ptcol="blue", ...)
{ 
  fhat <- x
  rm(x)
  
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

  ncont <- length(cont)
  
  if (missing(alphavec)) alphavec <- seq(0.1,0.5,length=ncont)
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

  if (!missing(y.group)) y.gr <- sort(unique(y.group))

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

  rgl.clear()
  rgl.bg(color="white")
  
  for (j in 1:m)
  {
    for (i in 1:ncont) 
      contour3d(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]],
                z=fhat$eval.points[[3]], f=fhat$estimate[[j]],
                level=hts[ncont-i+1],
                add=TRUE, alpha=alphavec[i], col=colors[j],...)

    if (drawpoints)   ## plot points
    {
      if (missing(y))
        points3d.rh(fhat$x[[j]][,1], fhat$x[[j]][,2], fhat$x[[j]][,3],
                    col=ptcol[j], size=size)
      else
      {
        if (missing(y.group))
          ppoints3d.rh(y[,1], y[,2], y[,3], col=ptcol[j], size=size)
        else
        {
          y.temp <- y[y.group==y.gr[j],]
          points3d.rh(y.temp[,1], y.temp[,2], y.temp[,3], col=ptcol[j], size=size)
        }
      }
    }
  }
  
  lines3d(c(origin[1],endpts[1]),rep(origin[2],2),rep(origin[3],2),size=3,
          color="black", add=TRUE)
  lines3d(rep(origin[1],2),c(origin[2],endpts[2]),rep(origin[3],2),size=3,
          color="black", add=TRUE)
  lines3d(rep(origin[1],2),rep(origin[2],2),c(origin[3],endpts[3]),size=3,
          color="black",add=TRUE)

  texts3d.rh(endpts[1]+0.1*abs(endpts[1]),origin[2],origin[3],xlabs,color="black",size=3)
  texts3d.rh(origin[1],endpts[2]+0.1*abs(endpts[2]),origin[3],ylabs,color="black",size=3)
  texts3d.rh(origin[1],origin[2],endpts[3]+0.1*abs(endpts[3]),zlabs,color="black",size=3)
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
  if (missing(xlim))
    xlim <- range(x[,1])
  if (missing(ylim))
    ylim <- range(x[,2])
  if (missing(zlim) & d==3)
    zlim <- range(x[,3])
  if (missing(gridsize)) 
    gridsize <- rep(100,d)

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

  prior.prob <- rep(0, length(gr))
  for (j in 1:length(gr))
    prior.prob[j] <- length(which(x.group==gr[j]))
  prior.prob <- prior.prob/nrow(x)
  fhat$prior.prob <- prior.prob
  
  fhat$type <- type
  
  class(fhat) <- "dade"
  
  return(fhat)
}

