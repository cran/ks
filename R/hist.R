#############################################################################
## Histogram density estimators
#############################################################################

histde <- function(x, binw, xmin, xmax, adj=0) 
{
    if (is.vector(x))
    {
        d <- 1; n <- length(x)
        if (missing(binw))
        binw <- 2*3^(1/(d+2))*pi^(2/(2*d+4))*sd(x)*length(x)^(-1/(d+2))
        nbin <- round(diff(range(x))*1.2/binw, 0)
    }
    else
    {
        d <- ncol(x); n <- nrow(x)
        if (missing(binw))
            binw <- 2*3^(1/(d+2))*pi^(2/(2*d+4))*apply(x, 2, sd)*nrow(x)^(-1/(d+2))
        nbin <- round(apply(apply(x, 2, range), 2, diff)*1.2/binw,0)     
    }
    
    if (d==1) fhat <- hist.1d(x, binw=binw, xmin=xmin, xmax=xmax, adj=adj)
    else if (d==2) fhat <- hist.2d(x=x, binw=binw, xmin=xmin, xmax=xmax, adj=adj)

    fhat$names <- parse.name(x)
    class(fhat) <- "histde"

    return(fhat)
 }

## 1D histogram

hist.1d <- function(x, nbin, binw, x.cut, xmin, xmax, adj=0, ...)
{
    if (missing(xmin)) xmin <- min(x)
    if (missing(xmax)) xmax <- max(x)
    if (missing(nbin)) nbin <- round(diff(range(x))*1.2/binw, 0)
    if (missing(x.cut)) x.cut <- seq(from=xmin-0.1*(xmax-xmin), to=xmax+0.1*(xmax-xmin), length=nbin+1)
    if (missing(binw)) binw <- diff(x.cut)[1]
    x.cut <- x.cut + adj*binw

    hs <- hist(x=x, breaks=x.cut, plot=FALSE) 
    hs <- list(x=x, estimate=hs$density, eval.points=hs$breaks, binw=binw, nbin=nbin)
    class(hs) <- "histde"

    return(hs)
}

## 2D histogram

hist.2d <- function(x, nbin, binw, x.cut, xmin, xmax, adj=0, ...)
{ 
    if (missing(nbin)) nbin <- round(apply(apply(x, 2, range), 2, diff)*1.2/binw,0)
    if (length(nbin)==1) nbin <- rep(nbin,2)
    if (missing(xmin)) xmin <- apply(x, 2, min)
    if (missing(xmax)) xmax <- apply(x, 2, max)
    if (missing(x.cut)) x.cut <- list(seq(from=xmin[1]-0.1*(xmax[1]-xmin[1]), to=xmax[1]+0.1*(xmax[1]-xmin[1]), length=nbin[1]+1), seq(from=xmin[2]-0.1*(xmax[2]-xmin[2]), to=xmax[2]+0.1*(xmax[2]-xmin[2]), length=nbin[2]+1))
    if (missing(binw)) binw <- c(diff(x.cut[[1]])[1], diff(x.cut[[2]])[1])
    x.cut[[1]] <- x.cut[[1]] + adj*binw[1]
    x.cut[[2]] <- x.cut[[2]] + adj*binw[2] 
    index.x <- cut(x[,1], x.cut[[1]], include.lowest=TRUE)
    index.y <- cut(x[,2], x.cut[[2]], include.lowest=TRUE)

    m <- matrix(0, nrow=nbin[1], ncol=nbin[2], dimnames=list(levels(index.x), levels(index.y)))
    for (i in 1:length(index.x))
        m[index.x[i], index.y[i]] <- m[index.x[i], index.y[i]] + 1

    hs2d <- list(x=x, estimate=m/(nrow(x)*prod(binw)), eval.points=x.cut, binw=binw, nbin=nbin)
    class(hs2d) <- "histde"

    return(hs2d)
}

#############################################################################
## plot histograms
#############################################################################

plot.histde <- function(x, ...)
{
    if (is.vector(x$x)) plot.histde.1d(fhat=x, ...)
    else plot.histde.2d(fhat=x, ...)
    invisible()
}

plot.histde.1d <- function(fhat, xlab, ylab="Density function", add=FALSE, drawpoints=FALSE, col="transparent", col.pt=4, jitter=FALSE, border=1, alpha=1, ...)
{
    if (missing(xlab)) xlab <- fhat$names
    col <- transparency.col(col, alpha=alpha)
    if (!add) plot(fhat$eval.points, c(fhat$estimate,0), type="n",  xlab=xlab, ylab=ylab, ...)
    rect(fhat$eval.points[-length(fhat$eval.points)], 0, fhat$eval.points[-1], fhat$estimate, border=border, col=col, ...)
    if (drawpoints)
        if (jitter) rug(jitter(fhat$x), col=col.pt)
        else rug(fhat$x, col=col.pt)

}

plot.histde.2d <- function(fhat, breaks, nbreaks=11, xlab, ylab, zlab="Density function", cex=1, pch=1, add=FALSE, drawpoints=FALSE, col, col.fun, alpha=1, col.pt=4, lty.rect=2, cex.text=1, border, lwd.rect=1, col.rect="transparent", add.grid=TRUE, ...)
{
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (missing(border)) border <- grey(0.5)
    if (!add) plot(fhat$x, col=col.pt, type="n", xlab=xlab, ylab=ylab, ...)
    if (missing(breaks)) breaks <- seq(min(fhat$estimate,0), max(fhat$estimate)+0.1*diff(range(fhat$estimate)), length=nbreaks)
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="heat", rev=TRUE, alpha=alpha) } 
    if (missing(col)) col <- col.fun(n=length(breaks))
    col <- transparency.col(col, alpha=alpha)
    
    for (i in 1:(nrow(fhat$estimate)))
        for (j in 1:(ncol(fhat$estimate)))
        {
            if (fhat$estimate[i,j]>=breaks[2])
            {
                rect(fhat$eval.points[[1]][i], fhat$eval.points[[2]][j], fhat$eval.points[[1]][i+1], fhat$eval.points[[2]][j+1], col=col[findInterval(fhat$estimate[i,j], breaks)], lty=lty.rect, border=border, lwd=lwd.rect)
            }
        }
    if (add.grid)
    {
        for (i in 1:length(fhat$eval.points[[1]]))
            lines(rep(fhat$eval.points[[1]][i],2), range(fhat$eval.points[[2]]), col=border, ...)
        for (j in 1:length(fhat$eval.points[[2]]))
            lines(range(fhat$eval.points[[1]]), rep(fhat$eval.points[[2]][j],2), col=border, ...)
    }
    
    if (drawpoints) 
        points(fhat$x[,1], fhat$x[,2], col=col.pt, cex=cex, pch=pch)
}

## predict method

predict.histde <- function(object, ..., x)
{
    fhat <- object$estimate
    if (!is.list(object$eval.points)) 
        d <- 1
    else d <- length(object$eval.points)
    if (d==1) 
    {
        gs <- length(object$eval.points)
        x.ind <- findInterval(x, object$eval.points, all.inside = FALSE)
        fhat[x.ind==0] <- object$estimate[1]
        fhat[x.ind==gs] <- object$estimate[gs]
    }
    else 
    {
        x <- matrix(x, ncol = d)
        gs <- sapply(object$eval.points, length)
        x.ind <- matrix(0, nrow=nrow(x), ncol=d)
        for (i in 1:d) x.ind[, i] <- findInterval(x[, i], object$eval.points[[i]], all.inside = FALSE)
        x.ind[x.ind==0] <- 1
        x.ind.flag <- x.ind==1
        for (i in 1:d) x.ind.flag[, i] <- x.ind.flag[, i] | x.ind[, i] ==gs[i]
    }

    return(fhat[x.ind])
}

## contourLevels method

contourLevels.histde <- function(...)
{
    return(contourLevels.kde(...))
}
