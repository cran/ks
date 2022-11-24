#############################################################################
## Kernel copula and copula density estimators
#############################################################################

## empirical pseudo-uniform transformation
## taken from pobs() function in copula package

pseudo.unif.empirical <- function (x, y) ##, na.last="keep", ties.method="average")
{
    if (missing(y)) y <- x
    if (is.vector(y)) y <- matrix(y, nrow=1)
    d <- ncol(x)
    u <- matrix(0, ncol=d, nrow=nrow(x))
    for (i in 1:d)
    {
        ecdf.fun <- ecdf(x=x[,i])
        u[,i] <- ecdf.fun(y[,i])
    }

    return(u)
}

## kernel pseudo-uniform transformation

pseudo.unif.kernel <- function(x, y, hs, binned=TRUE)
{
    if (missing(y)) y <- x
    if (is.vector(y)) y <- matrix(y, nrow=1)
    d <- ncol(x)
    u <- list()
    for (i in 1:d)
    {
        u[[i]] <- kcde(x=x[,i], h=hs[i], eval.points=y[,i], binned=binned)
    }
    u2 <- sapply(u, getElement, "estimate")

    return(u2)
}

#############################################################################
## Kernel copula estimator
#############################################################################

kcopula <- function(x, H, hs, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, marginal="kernel", verbose=FALSE)
{
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize

    ## clip data to xmin,xmax grid 
    grid.clip <- binned
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }
  
    ## default bandwidths
    if (missing(H)) H <- Hpi.kcde(x=x, binned=default.bflag(d=d,n=n))
    if (missing(hs))
    {
        hs <- rep(0, d)
        for (i in 1:d) hs[i] <- hpi.kcde(x=x[,i], binned=TRUE)
    }

    Fhat <- kcde(x=x, H=H, gridsize=gridsize, binned=binned, bgridsize=bgridsize, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, w=w, tail.flag="lower.tail", verbose=verbose)
    xlims <- sapply(Fhat$eval.points, range)
    xlims[1,] <- xlims[1,] - 0.1*abs(apply(xlims, 2, diff))
    xlims[2,] <- xlims[2,] + 0.1*abs(apply(xlims, 2, diff))

    ## generate pseudo-uniform values
    marginal <- match.arg(marginal, c("kernel", "empirical")) 
    if (marginal=="kernel")
    {  
        ## kernel pseudo-uniform
        u <- list()
        u.eval.points <- list()
        for (i in 1:d)
        {
            u.eval.points[[i]] <- kcde(x=x[,i], h=hs[i], eval.points=Fhat$eval.points[[i]], xmin=xlims[1,i], xmax=xlims[2,i], binned=TRUE)
        }
        y <- pseudo.unif.kernel(x=x, y=x, hs=hs, binned=TRUE)
        ep <- lapply(u.eval.points, getElement, "estimate")
    }
    else if (marginal=="empirical")
    {
        ## empirical pseudo-uniform
        y <- pseudo.unif.empirical(x=x, y=x)
        ep <- lapply(1:d, function(i) { f <- ecdf(x=x[,i]); f(Fhat$eval.points[[i]]) })
        ep <- numeric()
        for (i in 1:d)
        {
            f <- ecdf(x=x[,i])
            ep <- c(ep, list(f(Fhat$eval.points[[i]])))
        }
    }

    Chat <- Fhat
    Chat$x <- y 
    Chat$x.orig <- x
    Chat$eval.points <- ep
    Chat$hs <- hs

    ## loess smoothing on a uniform grid
    if (d==2)
    {
        ## select smaller grid for d==2 for memory usage reasons
        subselect <- round(cbind(seq(1,length(Chat$eval.points[[1]]), length=101), seq(1,length(Chat$eval.points[[2]]), length=101)),0)
        eval.points.df <- data.frame(expand.grid(Chat$eval.points[[1]][subselect[,1]], Chat$eval.points[[2]][subselect[,2]]))
        names(eval.points.df) <- paste("x", 1:ncol(eval.points.df), sep="")
        eval.points.df <- data.frame(estimate=as.numeric(Chat$estimate[subselect[,1], subselect[,2]]), eval.points.df)
    }
    else if (d==3)
    {
        ## select smaller grid for d==3 for memory usage reasons
        subselect <- round(cbind(seq(1,length(Chat$eval.points[[1]]), length=21), seq(1,length(Chat$eval.points[[2]]), length=21), seq(1,length(Chat$eval.points[[3]]), length=21)),0)
        eval.points.df <- data.frame(expand.grid(Chat$eval.points[[1]][subselect[,1]], Chat$eval.points[[2]][subselect[,2]], Chat$eval.points[[3]][subselect[,3]]))
        names(eval.points.df) <- paste("x", 1:ncol(eval.points.df), sep="")
        eval.points.df <- data.frame(estimate=as.numeric(Chat$estimate[subselect[,1], subselect[,2], subselect[,3]]), eval.points.df)
    }

    if (d==2) Chat.loess <- loess(estimate ~ x1+x2, data=eval.points.df, span=0.1)
    else if (d==3)  Chat.loess <- loess(estimate ~ x1+x2+x3, data=eval.points.df, span=0.1)

    u.eval.points.regular <- list()
    for (i in 1:d)
    u.eval.points.regular[[i]] <- seq(0,1,length=length(Chat$eval.points[[i]]))
    u.eval.points.regular.df <- data.frame(expand.grid(u.eval.points.regular))
    names(u.eval.points.regular.df) <- paste("x", 1:ncol(u.eval.points.regular.df), sep="")

    Chat.smoothed <- Chat
    Chat.smoothed$eval.points <- u.eval.points.regular
    Chat.smoothed$estimate <- array(predict(Chat.loess, newdata=u.eval.points.regular.df), dim=dim(Chat$estimate))

    ## interpolate NA boundary values from loess smoothing
    gsdim <- dim(Chat$estimate)
    if (d==2)
    {
        Chat.smoothed$estimate[1,] <- 0
        Chat.smoothed$estimate[gsdim[1],] <- Chat.smoothed$estimate[gsdim[1]-1,]*1.001
        Chat.smoothed$estimate[,1] <-0
        Chat.smoothed$estimate[,gsdim[2]] <- Chat.smoothed$estimate[,gsdim[2]-1]*1.001
        Chat.smoothed$estimate[gsdim[1],gsdim[2]] <- 1
        Chat.smoothed$estimate[Chat.smoothed$estimate>1] <- 1
    }
    else if (d==3)
    {
        Chat.smoothed$estimate[,,1] <- 0
        for (k in 2:(gsdim[3]-1))
        {
           Chat.smoothed$estimate[,,k][1,] <- 0
           Chat.smoothed$estimate[,,k][gsdim[1],] <- Chat.smoothed$estimate[,,k][gsdim[1]-1,]*1.001
           Chat.smoothed$estimate[,,k][,1] <-0
           Chat.smoothed$estimate[,,k][,gsdim[2]] <- Chat.smoothed$estimate[,,k][,gsdim[2]-1]*1.001
         }
        Chat.smoothed$estimate[,,gsdim[3]] <- Chat.smoothed$estimate[,,gsdim[3]-1]*1.001
        Chat.smoothed$estimate[gsdim[1],gsdim[2],gsdim[3]] <- 1
        Chat.smoothed$estimate[Chat.smoothed$estimate>1] <- 1
    }

    Chat <- Chat.smoothed
    Chat$marginal <- marginal
    class(Chat) <- "kcopula"

    return(Chat)
}

#############################################################################
## Kernel copula density estimator
#############################################################################

kcopula.de <- function(x, H, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, marginal="kernel", boundary.supp, boundary.kernel="beta", verbose=FALSE)
{
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    if (missing(binned)) binned <- ksd$binned
    if (missing(bgridsize)) bgridsize <- ksd$bgridsize
    if (missing(gridsize)) gridsize <- ksd$gridsize

    ## clip data to xmin,xmax grid for binned estimation
    grid.clip <- binned
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }

    ## default bandwidths
    hs <- rep(0, d)
    for (i in 1:d) hs[i] <- hpi.kcde(x=x[,i], binned=TRUE)

    ## generate pseudo-uniform values
    marginal <- match.arg(marginal, c("kernel", "empirical")) 
    if (marginal=="kernel") 
        y <- pseudo.unif.kernel(x=x, y=x, hs=hs, binned=TRUE)
    else if (marginal=="empirical") 
        y <- pseudo.unif.empirical(x=x, y=x)  

    colnames(y) <- colnames(x)
    if (missing(H)) H <- Hns(y)

    ## kernel copula density is boundary kernel estimator
    if (d==2 | d==3) 
        chat <- kde.boundary(x=y, H=H, gridsize=gridsize, supp=supp, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=FALSE, boundary.kernel=boundary.kernel, verbose=verbose)
    else
        stop("kcopula.de requires 2-d or 3-d data.")

    ## normalise KDE to integrate to 1 
    chat$estimate <- chat$estimate/sum(chat$estimate*apply(sapply(chat$eval.points, diff), 1, prod)[1])
    chat$names <- parse.name(x) 
    chat$x.orig <- x
    chat$hs <- hs

    ## compute prob contour levels
    if (compute.cont & missing(eval.points))
    chat$cont <- contourLevels(chat, cont=1:99, approx=approx.cont)
    chat$marginal <- marginal
    class(chat) <- "kcopula.de"

    return(chat)
}

#############################################################################
## S3 methods
#############################################################################

## plot methods

plot.kcopula <- function(x, ...)
{
    plot.kcde(x, ...)
}

plot.kcopula.de <- function(x, ...)
{
    plot.kde(x, ...)
}

## predict methods

predict.kcopula <- function(object, ..., x, u)
{
    if (missing(u)) { if (object$marginal=="kernel") u <- pseudo.unif.kernel(x=object$x.orig, y=x, hs=object$hs) }
    return(predict.kde(object, ..., x=u))
}

predict.kcopula.de <- function(object, ..., x, u)
{
    if (missing(u)) { if (object$marginal=="kernel") u <- pseudo.unif.kernel(x=object$x.orig, y=x, hs=object$hs) }
    return(predict.kde(object, ..., x=u))
}

## contourLevel method

contourLevels.kcopula.de <- function(x, ...)
{
    x1 <- x; class(x1) <- "kde"  
    return(contourLevels(x=x1, ...))
}
