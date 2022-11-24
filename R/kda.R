##############################################################################
## Kernel discriminant analysis
###############################################################################

###############################################################################
## KDEs of individual densities for KDA - 1- to 3-dim
##
## Parameters
## x - data values
## group - group variable
## Hs - bandwidth matrices
##
## Returns
## List with components (class dade)
## x - list of data values
## eval.points - evaluation points of dnesity estimate
## estimate - list of density estimate
## H - list of bandwidth matrices
##############################################################################

kda <- function(x, x.group, Hs, hs, prior.prob=NULL, gridsize, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, kde.flag=TRUE)
{
    if (missing(eval.points)) eval.points <- x
    if (!is.factor(x.group)) x.group <- factor(x.group) 
    gr <- levels(x.group)
    m <- length(gr)
    
    ## default values 
    ksd <- ks.defaults(x=x, w=w, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize  

    if (!is.null(prior.prob))
        if (!(identical(all.equal(sum(prior.prob), 1), TRUE))) stop("Sum of weights not equal to 1")
    
    ## clip data to xmin,xmax grid for binned estimation
    grid.clip <- binned 
    if (grid.clip)
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        xt <- truncate.grid(x=x, y=1:n, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- w[xt$y]; x.group <- x.group[xt$y]; n <- length(w)
    }

    if (d==1)
    {
        if (missing(hs)) hs <- hkda(x=x, x.group=x.group, bw="plugin", nstage=2, binned=default.bflag(d=d,n=n))
        
        ## compute KDA on grid
        if (kde.flag)
            fhat.list <- kda.1d(x=x, x.group=x.group, hs=hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=binned, bgridsize=bgridsize, xmin=xmin, xmax=xmax, compute.cont=compute.cont, approx.cont=approx.cont, w=w)
        
        ## compute KDA at eval.points
        fhat <- kda.1d(x=x, x.group=x.group, hs=hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=binned, bgridsize=bgridsize, xmin=xmin, xmax=xmax, eval.points=eval.points, compute.cont=compute.cont, approx.cont=approx.cont, w=w)
        fhat.wt <- matrix(0, ncol=m, nrow=length(eval.points)) 
    }
    else
    {
        if (missing(Hs)) Hs <- Hkda(x=x, x.group=x.group, bw="plugin", binned=default.bflag(d=d, n=n))
        
        ## Compute KDA on grid
        if (d>3) kde.flag <- FALSE
        if (kde.flag)
            fhat.list <- kda.nd(x=x, x.group=x.group, Hs=Hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=binned, bgridsize=bgridsize, xmin=xmin, xmax=xmax, compute.cont=compute.cont, approx.cont=approx.cont)
        
        ## compute KDA at eval.points
        fhat <- kda.nd(x=x, x.group=x.group, Hs=Hs, prior.prob=prior.prob, gridsize=gridsize, supp=supp, binned=FALSE, bgridsize=bgridsize, xmin=xmin, xmax=xmax, eval.points=eval.points, compute.cont=compute.cont, approx.cont=approx.cont)
        fhat.wt <- matrix(0, ncol=m, nrow=nrow(eval.points))  
    }

    for (j in 1:m)
        fhat.wt[,j] <- fhat$estimate[[j]]* fhat$prior.prob[j]
    
    ## assign y according largest weighted density value 
    disc.gr.temp <- apply(fhat.wt, 1, which.max)
    disc.gr <- as.factor(gr[disc.gr.temp])
    if (is.numeric(gr)) disc.gr <- as.numeric(levels(disc.gr))[disc.gr]
    
    if (kde.flag) fhat.list$x.group.estimate <- disc.gr
    else fhat.list <- disc.gr
    fhat.list$type <- "kda"
    
    return(fhat.list)
}

kda.1d <- function(x, x.group, hs, prior.prob, gridsize, supp, eval.points, binned, bgridsize, xmin, xmax, w, compute.cont, approx.cont)
{
    gr <- levels(x.group)
    m <- length(gr)
    d <- 1
    hmax <- max(hs)
    if (missing(xmin)) xmin <- min(x) - supp*hmax
    if (missing(xmax)) xmax <- max(x) + supp*hmax
    fhat.list <- list()
    
    for (j in 1:m)
    {
        xx <- x[x.group==gr[j]]
        ww <- w[x.group==gr[j]]
        h <- hs[j]
        
        ## compute individual density estimate
        if (missing(eval.points))
            fhat.temp <- kde(x=xx, h=h, supp=supp, xmin=xmin, xmax=xmax, bgridsize=bgridsize, gridsize=gridsize, w=ww, binned=binned)
        else
            fhat.temp <- kde(x=xx, h=h, w=ww, binned=binned, bgridsize=bgridsize, gridsize=gridsize, eval.points=eval.points)
    
        fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$estimate))
        fhat.list$eval.points <- fhat.temp$eval.points
        fhat.list$x <- c(fhat.list$x, list(xx))
        fhat.list$h <- c(fhat.list$h, h)
        fhat.list$H <- c(fhat.list$H, h^2)
        fhat.list$w <- c(fhat.list$w, list(ww))
        
        ## compute prob contour levels
        if (compute.cont & missing(eval.points))
        {
            contlev <- contourLevels(fhat.temp, cont=1:99, approx.cont=approx.cont)
            fhat.list$cont <- c(fhat.list$cont, list(contlev))
        }
    }
    fhat.list$names <- parse.name(x)  ## add variable names
    fhat.list$binned <- binned
    fhat.list$gridded <- fhat.temp$gridded
    
    if (is.null(prior.prob))
    {
        pr <- rep(0, length(gr))
        for (j in 1:length(gr)) pr[j] <- length(which(x.group==gr[j]))
        pr <- pr/length(x)
        fhat.list$prior.prob <- pr
    }
    else
        fhat.list$prior.prob <- prior.prob
    
    fhat.list$x.group <- x.group
    class(fhat.list) <- "kda"
    
    return(fhat.list) 
}

kda.nd <- function(x, x.group, Hs, prior.prob, gridsize, supp, eval.points, binned, bgridsize, xmin, xmax, w, compute.cont, approx.cont)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    gr <- levels(x.group)
    m <- length(gr)
    d <- ncol(x)
  
    ## find largest bandwidth matrix to initialise grid
    detH <- vector() 
    for (j in 1:m)
        detH[j] <- det(Hs[((j-1)*d+1) : (j*d),])  
    Hmax.ind <- which.max(detH)
    Hmax <- Hs[((Hmax.ind-1)*d+1) : (Hmax.ind*d),]
    if (missing(xmin)) xmin <- apply(x, 2, min) - supp*max(sqrt(diag(Hmax)))
    if (missing(xmax)) xmax <- apply(x, 2, max) + supp*max(sqrt(diag(Hmax)))
    if (missing(w)) w <- rep(1, nrow(x))

    if (binned & d > 4) stop("Binning only available for 1- to 4-d data")
    if (missing(bgridsize)) bgridsize <- default.gridsize(d)
    if (missing(gridsize)) gridsize <- default.gridsize(d)

    fhat.list <- list()
    for (j in 1:m)
    {
        xx <- x[x.group==gr[j],]
        ww <- w[x.group==gr[j]]   
        H <- Hs[((j-1)*d+1) : (j*d),]     
        
        ## compute individual density estimate
        if (binned)
            fhat.temp <- kdde(x=xx, bgridsize=bgridsize, H=H, xmin=xmin, xmax=xmax, w=ww, deriv.order=0, binned=TRUE)
        else if (missing(eval.points))
            fhat.temp <- kde(x=xx, H=H, supp=supp, xmin=xmin, xmax=xmax, gridsize=gridsize, w=ww)
        else
            fhat.temp <- kde(x=xx, H=H, eval.points=eval.points, w=ww)
        
        fhat.list$estimate <- c(fhat.list$estimate, list(fhat.temp$estimate))
        fhat.list$eval.points <- fhat.temp$eval.points
        fhat.list$x <- c(fhat.list$x, list(xx))
        fhat.list$H <- c(fhat.list$H, list(H))
        fhat.list$w <- c(fhat.list$w, list(ww))
        
        ## compute prob contour levels
        if (compute.cont & missing(eval.points))
        {
            contlev <- contourLevels(fhat.temp, cont=1:99, approx.cont=approx.cont)
            fhat.list$cont <- c(fhat.list$cont, list(contlev))
        }
    }

    fhat.list$names <- parse.name(x)  ## add variable names
    fhat.list$binned <- binned
    fhat.list$gridded <- fhat.temp$gridded

    if (is.null(prior.prob))
    {
        pr <- rep(0, length(gr))
        for (j in 1:length(gr)) pr[j] <- length(which(x.group==gr[j]))
        pr <- pr/nrow(x)
        fhat.list$prior.prob <- pr
    }
    else
        fhat.list$prior.prob <- prior.prob

    fhat.list$x.group <- x.group
    class(fhat.list) <- "kda"

    return (fhat.list)
}

  
##############################################################################
## S3 methods for kda objects
##############################################################################

## contourLevel method

contourLevels.kda <- function(x, prob, cont, nlevels=5, approx=TRUE,...) 
{
    fhat <- x
    m <- length(fhat$x)
    hts <- list()

    for (j in 1:m)
    {
        fhatj <- list(x=fhat$x[[j]], eval.points=fhat$eval.points, estimate=fhat$estimate[[j]], H=fhat$H[[j]], binned=fhat$binned, gridded=fhat$gridded)
        class(fhatj) <- "kde"
        hts[[j]] <- contourLevels(x=fhatj, prob=prob, cont=cont, nlevels=nlevels, approx=approx, ...)
    }
   
    return(hts) 
}

## predict method

predict.kda <- function(object, ..., x)
{
    fhat <- object
    m <- length(fhat$prior.prob)
    if (is.vector(fhat$x[[1]])) n <- length(x) else { if (is.vector(x)) n <- 1 else n <- nrow(x) }
    fhat.temp <- matrix(0, ncol=m, nrow=n)
    for (j in 1:m)    
        fhat.temp[,j] <- fhat$prior.prob[j]*grid.interp(x=x, gridx=fhat$eval.points, f=fhat$estimate[[j]])
    
    est.group <- apply(fhat.temp, 1, which.max)
    est.group <- as.factor(sort(unique(fhat$x.group))[est.group])
    
    return(est.group)
}

## plot method
## fhat - output from `kda.kde'
## y - data points (separate from training data inside fhat)
## y.group - data group labels

plot.kda <- function(x, y, y.group, ...) 
{
    if (is.vector(x$x[[1]]))
        plotkda.1d(x=x, y=y, y.group=y.group,  ...)
    else
    {  
        d <- ncol(x$x[[1]])      
        if (d==2)
        {
            opr <- options()$preferRaster; if (!is.null(opr)) if (!opr) options("preferRaster"=TRUE)
            plotkda.2d(x=x, y=y, y.group=y.group, ...)
            if (!is.null(opr)) options("preferRaster"=opr) 
        }
        else if (d==3)  
            plotkda.3d(x=x, y=y, y.group=y.group, ...)
    }
}

plotkda.1d <- function(x, y, y.group, prior.prob=NULL, xlim, ylim, xlab, ylab="Weighted density function", drawpoints=FALSE, col, col.fun, col.part, col.pt, lty, jitter=TRUE, rugsize, add=FALSE, alpha=1, ...)
{ 
    fhat <- x
    m <- length(fhat$x)
    if (missing(xlab)) xlab <- fhat$names
    if (is.null(prior.prob))
    prior.prob <- fhat$prior.prob

    if (m != length(prior.prob))
        stop("prior.prob not same length as number of components in fhat")
    if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
        stop("Sum of weights not equal to 1")

    weighted.fhat <- matrix(0, nrow=length(fhat$eval.points), ncol=m) 
    for (j in 1:m)  weighted.fhat[,j] <- fhat$estimate[[j]]*fhat$prior.prob[j]

    if (missing(xlim)) xlim <- range(fhat$eval.points)
    if (missing(ylim)) ylim <- range(weighted.fhat)
    if (missing(lty)) lty <- rep(1, m)
    if (length(lty) < m) lty <- rep(lty, m)
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Dark2") }
    if (missing(col)) col <- col.fun(m)
    if (length(col) < m) col <- rep(col, m)
    if (missing(col.part)) col.part <- plot3D::alpha.col(col, alpha=alpha) 
    if (missing(col.pt)) col.pt <- col.fun(m)
    if (length(col.pt)==1) col.pt <- rep(col.pt, m)
  
    ## plot each training group's KDE in separate colour and line type 
    if (!add) plot(fhat$eval.points, weighted.fhat[,1], type="l", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, lty=lty[1], col=col[1], ...)
  
    if (m > 1)
        for (j in 2:m) lines(fhat$eval.points, weighted.fhat[,j], lty=lty[j], col=col[j], ...)

    ydata <- seq(min(fhat$eval.points), max(fhat$eval.points), length=401)
    x.gr <- 1:m
    ydata.gr <- x.gr[apply(weighted.fhat,1, which.max)] 

    ## draw partition class as rug-like plot
    plot.lim <- par()$usr
    if (missing(rugsize)) rugsize <- abs(plot.lim[4]-plot.lim[3])/50

    image(ydata, c(plot.lim[3], plot.lim[3]+rugsize), cbind(as.numeric(ydata.gr), as.numeric(ydata.gr)), breaks=0.5+(0:m), add=TRUE, col=col.part, ...)  

    for (j in 1:m)
    {     
        ## draw data points
        if (drawpoints)
        {
            if (missing(y)) yy <- fhat$x[[j]]
            else 
            {
                if (missing(y.group)) yy <- fhat$x[[j]]
                else yy <- y[y.group==levels(y.group)[j]]
            }
            if (jitter) yy <- jitter(yy)
            rug(yy, col=col.pt[j], ticksize=-0.03)
        }
    }   
}

plotkda.2d <- function(x, y, y.group, prior.prob=NULL, display.part="filled.contour",
    cont=c(25,50,75), abs.cont, approx.cont=TRUE, xlim, ylim, xlab, ylab,
    drawpoints=FALSE, drawlabels=TRUE, cex=1, pch, lty, part=TRUE, col, col.fun, col.part, col.pt, alpha=1, lwd=1, lwd.part=0, labcex=1, add=FALSE, ...)
{ 
    fhat <- x
    m <- length(fhat$x)
    xtemp <- numeric()
    for (j in 1:m) xtemp <- rbind(xtemp, fhat$x[[j]]) 
    if (missing(xlim)) xlim <- range(xtemp[,1])
    if (missing(ylim)) ylim <- range(xtemp[,2])
    if (missing(pch)) pch <- 1:m
    if (missing(lty)) lty <- rep(1, m)
    if (length(lty) < m) lty <- rep(lty, m)
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Dark2", alpha=1) }
    if (missing(col)) col <- col.fun(m)
    if (length(col) < m) col <- rep(col, m)
    if (missing(col.part)) col.part <- plot3D::alpha.col(col, alpha=alpha) 
    if (missing(col.pt)) col.pt <- col.fun(m)
    if (length(col.pt)==1) col.pt <- rep(col.pt, m)

    x.names <- fhat$names
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
    if (is.null(prior.prob)) prior.prob <- fhat$prior.prob

    if (m != length(prior.prob))
        stop("prior.prob not same length as number of components in fhat")
    if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
        stop("Sum of weights not equal to 1")

    ## set up plot
    if (!add)
    {
        if (missing(y)) 
            plot(fhat$x[[1]], type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
        else
        plot(y, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
    }

    ## set up common grid for all densities 
    class.grid <- array(0, dim=dim(fhat$estimate[[1]]))
    temp <- matrix(0, ncol=length(fhat$est), nrow=nrow(fhat$est[[1]]))
    for (j in 1:ncol(fhat$estimate[[1]]))
    {
        for (k in 1:length(fhat$estimate))
            temp[,k] <- fhat$estimate[[k]][,j]* prior.prob[k]
        class.grid[,j] <- max.col(temp)
    }

    ## draw partition
    fhat.part <- fhat
    fhat.part$estimate <- class.grid  
    fhat.part$H <- fhat$H[[1]]
    fhat.part$x <- fhat$x[[1]]
    fhat.part$w <- fhat$w[[1]]
    fhat.part$cont <- fhat$cont[[1]]
    class(fhat.part) <- "kde.part"

    if (part) plot(fhat.part, col=col.part, add=TRUE, display=display.part, lwd=lwd.part, alpha=alpha, ...)

    ## common contour levels removed from >= v1.5.3 
    if (missing(abs.cont))
    {
        hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        nhts <- length(hts[[1]])
    }
    else
    {
        hts <- abs.cont
        nhts <- length(hts)
    }
 
        ## draw contours
        for (j in 1:m)
        {
            for (i in 1:nhts) 
            {
                if (missing(abs.cont))
                {
                    scale <- cont[i]/hts[[j]][i]
                    contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                        fhat$estimate[[j]]*scale, level=hts[[j]][i]*scale, add=TRUE, 
                        drawlabels=drawlabels, lty=lty[j], col=col[j], lwd=lwd, labcex=labcex, 
                        ...)
                }
                else
                {
                    contour(fhat$eval.points[[1]], fhat$eval.points[[2]], 
                        fhat$estimate[[j]], level=hts[i], add=TRUE, 
                        drawlabels=drawlabels, lty=lty[j], col=col[j], lwd=lwd, labcex=labcex, 
                        ...)
                }
            }
        }
  
    for (j in 1:m)
    {  
        ## draw data points
        if (drawpoints)
        {
            if (missing(y))
                points(fhat$x[[j]], pch=pch[j], col=col.pt[j], cex=cex)
            else 
            {
                if (missing(y.group))
                {  
                    points(y, col=col.pt[1], cex=cex)
                    j <- m+1
                }
                else
                    points(y[y.group==levels(y.group)[j],], pch=pch[j], col=col.pt[j], cex=cex) 
            }
        }
    }   
}

plotkda.3d <- function(x, y, y.group, prior.prob=NULL, display="plot3D", cont=c(25,50,75), abs.cont, approx.cont=TRUE, colors, col, col.fun, col.pt, alpha=0.5, alphavec, xlab, ylab, zlab, drawpoints=FALSE, size=3, cex=1, pch, theta=-30, phi=40, d=4, ticktype="detailed", bty="f", add=FALSE, ...)
{
    fhat <- x
    
    m <- length(fhat$x) 
    if (is.null(prior.prob)) prior.prob <- fhat$prior.prob
    if (m != length(prior.prob))
        stop("prior.prob not same length as number of components in fhat")
    if (!(identical(all.equal(sum(prior.prob), 1), TRUE)))  
        stop("Sum of prior weights not equal to 1")
    
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Dark2") }
    if (missing(col)) col <- col.fun(m)
    if (length(col) < m) col <- rep(col, m)
    if (missing(col.pt)) col.pt <- col.fun(m)
    if (length(col.pt)==1) col.pt <- rep(col.pt, m)
    if (missing(pch)) pch <- 1:m
    if (length(pch)==1) pch <- rep(pch, m)
   
    x.names <- colnames(fhat$x[[1]])
    if (missing(xlab)) if (is.null(x.names)) xlab <- "x" else xlab <- x.names[1]
    if (missing(ylab)) if (is.null(x.names)) ylab <- "y" else ylab <- x.names[2]
    if (missing(zlab)) if (is.null(x.names)) zlab <- "z" else zlab <- x.names[3]
    
    xx <- numeric(0)
    for (j in 1:m) xx <- rbind(xx, fhat$x[[j]])
    
    ## common contour levels removed from >= v1.5.3 
    if (missing(abs.cont))
    {
        hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        nhts <- length(hts[[1]])
    }
    else
    {
        hts <- abs.cont
        nhts <- length(hts)
    }
    
    if (missing(alphavec)) alphavec <- seq(0.05,0.2,length=nhts)
    
    disp <- match.arg(display, c("plot3D", "rgl")) 
    if (disp %in% "plot3D")
    {
        for (j in 1:m)
        {
            for (i in 1:nhts)
            { 
                cti <- hts[[j]][nhts-i+1]
                if (cti <= max(fhat$estimate[[j]]))
                    plot3D::isosurf3D(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], colvar=fhat$estimate[[j]], level=cti, add=add | (j>1 | i>1), alpha=alphavec[i], col=col[j], phi=phi, theta=theta, xlab=xlab, ylab=ylab, zlab=zlab, d=d, ticktype=ticktype, bty=bty, ...)
            }
        }
        
        ## plot points
        if (drawpoints)   
        {
            for (j in 1:m)
            {
                if (missing(y))
                    plot3D::points3D(fhat$x[[j]][,1], fhat$x[[j]][,2], fhat$x[[j]][,3], col=col.pt[j], cex=cex, alpha=1, add=TRUE, pch=pch[j])
                else
                {
                    if (missing(y.group))
                        plot3D::points3D(y[,1], y[,2], y[,3], col=col.pt, cex=cex, alpha=alpha, add=TRUE, pch=pch[j])
                    else
                    {
                        y.temp <- y[y.group==levels(y.group)[j],]
                        if (nrow(y.temp)>0)
                            plot3D::points3D(y.temp[,1], y.temp[,2], y.temp[,3], col=col.pt[j], cex=cex, alpha=alpha, add=TRUE, pch=pch[j])
                    }
                }
            }
        }
    
    }
    else if (disp %in% "rgl")
    {
        ## suggestions from Viktor Petukhov 08/03/2018
        if (!requireNamespace("rgl", quietly=TRUE)) stop("Install the rgl package as it is required.", call.=FALSE)
        if (!requireNamespace("misc3d", quietly=TRUE)) stop("Install the misc3d package as it is required.", call.=FALSE)

        if (missing(colors)) colors <- col
        
        xtemp <- numeric(); for (i in 1:length(fhat$x)) xtemp <- rbind(xtemp, fhat$x[[i]])
        rgl::plot3d(x=xtemp[,1], y=xtemp[,2], z=xtemp[,3], type="n", xlab=xlab, ylab=ylab, zlab=zlab, ...)

        for (j in 1:m)
        {
            for (i in 1:nhts)
            { 
                cti <- hts[[j]][nhts-i+1]
                if (cti <= max(fhat$estimate[[j]]))
                    misc3d::contour3d(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], f=fhat$estimate[[j]], level=cti, add=TRUE, alpha=alphavec[i], color=colors[j], ...)
            }
            if (drawpoints)   ## plot points
            {
                if (missing(y))
                    rgl::points3d(fhat$x[[j]][,1], fhat$x[[j]][,2], fhat$x[[j]][,3], color=col.pt[j], size=size, alpha=alpha)
                else
                {
                    if (missing(y.group))
                        rgl::points3d(y[,1], y[,2], y[,3], color=col.pt, size=size, alpha=1)
                    else
                    {
                        y.temp <- y[y.group==levels(y.group)[j],]
                        if (nrow(y.temp)>0)
                            rgl::points3d(y.temp[,1], y.temp[,2], y.temp[,3], color=col.pt[j], size=size, alpha=alpha)
                    }
                }
            }
        }
    }   
}

###############################################################################
## Find bandwidths for each class in training set, for 2- to 6-dim 
##
## Parameters
## x - data values
## group - group variable
## bw - type of bandwidth selector
## nstage, pilot, pre - parameters for plugin bandwidths
## diag - FALSE - use full b/w matrices
##      - TRUE - use diag b/w matrices
##
## Returns
## Matrix of bandwidths for each group in training set
###############################################################################

hkda <- function(x, x.group, bw="plugin", ...)
{
    gr <- sort(unique(x.group))
    m <- length(gr)
    bw <- match.arg(bw, c("lscv", "plugin", "scv"))
    hs <- numeric(0) 
    for (i in 1:m)
    {
        y <- x[x.group==gr[i]]
        if (bw=="plugin") h <- hpi(y, ...)
        else if ((bw=="lscv") | (bw=="ucv")) h <- hlscv(y, ...)
        else if (bw=="scv") h <- hscv(y, ...)
        hs <- c(hs, h)
    }

    return(hs)
}
   
Hkda <- function(x, x.group, Hstart, bw="plugin", ...)
{
    d <- ncol(x)
    gr <- sort(unique(x.group))
    m <- length(gr)
    bw <- match.arg(bw, c("lscv", "plugin", "scv", "ucv")) 
    Hs <- numeric(0)

    for (i in 1:m)
    {
        y <- x[x.group==gr[i],]
        if (!missing(Hstart)) 
        {
            Hstarty <- Hstart[((i-1)*d+1) : (i*d),]
            if ((bw=="lscv") | (bw=="ucv")) H <- Hlscv(y, Hstart=Hstarty, ...)
            else if (bw=="scv") H <- Hscv(y, Hstart=Hstarty, ...)
            else if (bw=="plugin") H <- Hpi(y, Hstart=Hstarty, ...)
        }
        else
        {
            if ((bw=="lscv") | (bw=="ucv")) H <- Hlscv(y, ...)
            else if (bw=="scv") H <- Hscv(y, ...)
            else if (bw=="plugin") H <- Hpi(y, ...)
        }
        Hs <- rbind(Hs, H)
    }

    return(Hs)   
}

Hkda.diag <- function(x, x.group, bw="plugin", ...)
{
    d <- ncol(x)
    gr <- sort(unique(x.group))
    m <- length(gr)
    bw <- match.arg(bw, c("lscv", "plugin", "scv", "ucv"))
    Hs <- numeric(0)

    for (i in 1:m)
    {
        y <- x[x.group==gr[i],]
        if ((bw=="lscv") | (bw=="ucv")) H <- Hlscv.diag(y, ...)
        else if (bw=="plugin") H <- Hpi.diag(y, ...)
        else if (bw=="scv") H <- Hscv.diag(y, ...)
        Hs <- rbind(Hs, H)
    }

    return(Hs)   
}

###############################################################################
## Compares true group classification with an estimated one
##
## Parameters
## group - true group variable
## est.group - estimated group variable
##
## Returns
## List with components
## comp - cross-classification table of groupings - true groups are the rows,
##        estimated groups are the columns
## error - total mis-classification rate
###############################################################################

compare <- function(x.group, est.group, by.group=FALSE)
{
    if (length(x.group)!=length(est.group))
    stop("Group label vectors not the same length")

    if (!is.factor(x.group)) x.group <- factor(x.group)
    grlab <- levels(x.group)
    m <- length(grlab)
    comp <- table(x.group, est.group)

    if (by.group)
    {
        er <- vector()
        for (i in 1:m)
          er[i] <- 1-comp[i,i]/rowSums(comp)[i]
        er <- matrix(er, ncol=1)
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

    if (nrow(comp)==2 & nrow(comp)==2)
    {
        TN <- comp[1,1]; FP <- comp[1,2]; FN <- comp[2,1]; TP <- comp[2,2]
        spec <- 1-(FP/(FP+TN))
        sens <- 1-(FN/(TP+FN))
        comp <- list(cross=comp, error=er, TP=TP, FP=FP, FN=FN, TN=TN, spec=spec, sens=sens) 
    }
    else comp <- list(cross=comp, error=er)

    return(comp)
}

###############################################################################
## Computes cross-validated misclassification rates (for use when test data is
## not independent of training data) for KDA
##
## Parameters
## x - training data
## x.group - group variable for x
## y - data values to be classified
## Hs - bandwidth matrices
## prior.prob - prior probabilities
##
## Returns
## List with components
## comp - cross-classification table of groupings - true groups are the rows,
##        estimated groups are the columns
## error - total mis-classification rate
###############################################################################

compare.kda.cv <- function(x, x.group, bw="plugin", prior.prob=NULL, Hstart, by.group=FALSE, verbose=FALSE, recompute=FALSE, ...)
{
    if (verbose) pb <- txtProgressBar()
    bw <- match.arg(bw, c("lscv", "plugin", "scv", "ucv"))

    ## 1-d
    if (is.vector(x))
    {
        n <- length(x)
        h <- hkda(x, x.group, bw=bw, ...)
        gr <- sort(unique(x.group)) 
        kda.cv.gr <- x.group

        for (i in 1:n)
        {
            h.mod <- h
            ## find group that x[i] belongs to 
            ind <- which(x.group[i]==gr)
            indx <- x.group==gr[ind]
            indx[i] <- FALSE

            if ((bw=="lscv") | (bw=="ucv")) h.temp <- hlscv(x[indx], , ...)
            else if (bw=="plugin") h.temp <- hpi(x[indx], , ...)
            else if (bw=="scv") h.temp <- hscv(x[indx], , ...)
            h.mod[ind] <- h.temp

            ## recompute KDA estimate of groups with x[i] excluded
            if (verbose) setTxtProgressBar(pb, i/n)
            kda.cv.gr[i] <- kda(x[-i], x.group[-i], hs=h.mod, eval.points=x, prior.prob=prior.prob, kde.flag=FALSE)[i]
        }
        if (verbose) close(pb)
        
        return(compare(x.group, kda.cv.gr, by.group=by.group)) 
    }

    ## multi-dimensional   
    n <- nrow(x)
    d <- ncol(x)

    if (!missing(Hstart))
        H <- Hkda(x, x.group, bw=bw, Hstart=Hstart, ...)
    else
        H <- Hkda(x, x.group, bw=bw, ...)

    ## classify data x using KDA rules based on x itself
    ## kda.group <- kda(x, x.group, Hs=H, y=x, prior.prob=prior.prob)
    ## comp <- compare(x.group, kda.group)

    gr <- sort(unique(x.group)) 
    kda.cv.gr <- x.group

    for (i in 1:n)
    {
        H.mod <- H
        ## find group that x[i] belongs to 
        ind <- which(x.group[i]==gr)
        indx <- x.group==gr[ind]
        indx[i] <- FALSE

        if (recompute)
        {
            ## compute b/w matrix for that group with x[i] excluded
            if (!missing(Hstart))
            {  
                Hstart.temp <- Hstart[((ind-1)*d+1):(ind*d),]

                if (bw=="plugin") H.temp <- Hpi(x[indx,], Hstart=Hstart.temp,  ...)
                else if (bw=="scv") H.temp <- Hscv(x[indx,],  Hstart=Hstart.temp, ...)
                else if ((bw=="lscv") | (bw=="ucv")) H.temp <- Hlscv(x[indx,],  Hstart=Hstart.temp, ...)
            }
            else
            {
                if (bw=="plugin") H.temp <- Hpi(x[indx,], ...)
                else if (bw=="scv") H.temp <- Hscv(x[indx,], ...)
                else if ((bw=="lscv") | (bw=="ucv")) H.temp <- Hlscv(x[indx,], ...) 
            }

            H.mod[((ind-1)*d+1):(ind*d),] <- H.temp
        }
        ## recompute KDA estimate of groups with x[i] excluded
        if (verbose) setTxtProgressBar(pb, i/n)  

        kda.cv.gr[i] <- kda(x[-i,], x.group[-i], Hs=H.mod, eval.points=x, prior.prob=prior.prob, kde.flag=FALSE)[i]
    }
    if (verbose) close(pb)

    return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}

###############################################################################
## Same as compare.kda.cv except uses diagonal b/w matrices
###############################################################################

compare.kda.diag.cv <- function(x, x.group, bw="plugin", prior.prob=NULL, by.group=FALSE, verbose=FALSE, recompute=FALSE, ...)
{
    if (is.vector(x))  return(compare.kda.cv(x=x, x.group=x.group, by.group=by.group, verbose=verbose, prior.prob=prior.prob, recompute=recompute, ...))
    n <- nrow(x)
    d <- ncol(x)

    H <- Hkda.diag(x, x.group, bw=bw, ...)

    if (verbose) pb <- txtProgressBar()
    bw <- match.arg(bw, c("lscv", "plugin", "scv", "ucv"))
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
            if (bw=="plugin")
                H.temp <- Hpi.diag(x[indx,], ...)
            else if ((bw=="lscv") | (bw=="ucv"))
                H.temp <- Hlscv.diag(x[indx,], ...)
            else if (bw=="scv")
                H.temp <- Hscv.diag(x[indx,], ...)
            H.mod[((ind-1)*d+1):(ind*d),] <- H.temp
        }
    
        if (verbose) setTxtProgressBar(pb, i/n) 
        kda.cv.gr[i] <- kda(x[-i,], x.group[-i], Hs=H.mod, eval.points=x, prior.prob=prior.prob, kde.flag=FALSE)[i]  
    }
    if (verbose) close(pb)
  
    return(compare(x.group, kda.cv.gr, by.group=by.group)) 
}
