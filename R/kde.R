###############################################################################
## Multivariate kernel density estimators
###############################################################################

###############################################################################
## Multivariate kernel density estimate using normal kernels
##
## Parameters
## x - points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
## eval.points - compute density estimate at these points (if missing
##            and dim(x) = 2, 3 compute density estimate over grid)  
## eval.levels - compute 3-D in 2-D slices taken at these level curves   
##
## Returns
## list with first d components with the points that the density
## estimate is evaluated at, and values of the density estimate 
##############################################################################

kde <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, positive=FALSE, adj.positive, w, compute.cont=TRUE, approx.cont=TRUE, unit.interval=FALSE, density=FALSE, verbose=FALSE)
{
    ## default values
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize
    if (binned & d>4) stop("Binned estimation for d>4 not implemented. Set binned=FALSE for exact estimation.")

    ## clip data to xmin, xmax grid limits for binned estimation
    grid.clip <- binned    
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        if (positive & missing(xmin)) { xmin <- rep(0,d) } 
        if (unit.interval) { if (missing(xmin)) xmin <- rep(0,d); if (missing(xmax)) xmax <- rep(1,d) }
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }
    
    ## default bandwidths
    if (d==1 & missing(h)) 
    {
        if (positive) x1 <- log(x) else x1 <- x
        if (unit.interval) x1 <- qnorm(x)
        h <- hpi(x=x1, nstage=2, binned=default.bflag(d=d, n=n), deriv.order=0)
    }
    if (d>1 & missing(H))
    {
        if (positive) x1 <- log(x) else x1 <- x
        H <- Hpi(x=x1, nstage=2, binned=default.bflag(d=d, n=n), deriv.order=0)
    }
    
    ## compute binned estimator
    if (binned)
    {
        if (positive)
        {
            if (d==1)
            {
                fhat <- kde.positive.1d(x=x, h=h, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, binned=binned, adj.positive=adj.positive)
            }
            else if (d==2)
            {
                fhat <- kde.positive.2d(x=x, H=H, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, binned=binned, adj.positive=adj.positive)
            }
        }
        else if (unit.interval)
        {
             fhat <- kde.unit.interval.1d(x=x, binned=binned, h=h)
        }
        else 
        {
            fhat <- kdde.binned(x=x, H=H, h=h, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.order=0, verbose=verbose)
        }
        
        if (!missing(eval.points))
        {
            fhat$estimate <- predict(fhat, x=eval.points)
            fhat$eval.points <- eval.points
        }
    }
    else
    {
        ## compute exact (non-binned) estimator   
        ## 1-dimensional    
        if (d==1)
        {
            if (missing(eval.points))
            {
                if (unit.interval)
                    fhat <- kde.unit.interval.1d(x=x, h=h, binned=FALSE)
                else if (positive)
                    fhat <- kde.positive.1d(x=x, h=h, xmin=xmin, xmax=xmax, w=w, binned=FALSE, adj.positive=adj.positive)
                else
                    fhat <- kde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, xmin=xmin, xmax=xmax, adj.positive=adj.positive, gridtype=gridtype, w=w)
            }
            else
                fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive, w=w)
        }
        ## multi-dimensional
        else
        {  
            if (is.data.frame(x)) x <- as.matrix(x)
            
            if (missing(eval.points))
            {
                if (d==2)
                {
                    if (positive)
                        fhat <- kde.positive.2d(x=x, H=H, gridsize=gridsize, xmin=xmin, xmax=xmax, w=w, binned=binned, adj.positive=adj.positive)
                    else
                        fhat <- kde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
                }
                else if (d==3)
                    fhat <- kde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose) 
                else 
                    fhat <- kde.grid.nd(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
            }
            else 
                fhat <- kde.points(x=x, H=H, eval.points=eval.points, w=w, verbose=verbose)     
        }
    }
    
    if (density) fhat$estimate[fhat$estimate<0] <- 0
    fhat$binned <- binned
    fhat$names <- parse.name(x)  ## add variable names
    fhat$w <- w
    fhat$type <- "kde" 
    class(fhat) <- "kde"
  
    ## compute prob contour levels
    if (compute.cont & missing(eval.points))
        fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)
    
    return(fhat)
}

###############################################################################
## Univariate kernel density estimate on a grid
###############################################################################

kde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype, w)
{
    if (missing(xmin)) xmin <- min(x) - h*supp
    if (missing(xmax)) xmax <- max(x) + h*supp
    if (missing(gridtype)) gridtype <- "linear"

    if (positive)
    {
        if (missing(adj.positive)) adj.positive <- abs(min(x))
        y <- log(x + adj.positive)  ## transform positive data x to real line

        gridx <- seq(max(0, xmin), xmax, length=gridsize)
        gridy <- log(gridx + adj.positive)
        gridtype.vec <- "linear" 
    }
    else
    {
        y <- x
        gridtype1 <- match.arg(gridtype, c("linear", "sqrt", "quantile", "exp")) 
        if (gridtype1=="linear")
        {
            gridy <- seq(xmin, xmax, length=gridsize)
        }
        else if (gridtype1=="sqrt")
        {
            gridy.temp <- seq(sign(xmin)*sqrt(abs(xmin)), sign(xmax)*sqrt(abs(xmax)), length=gridsize)
            gridy <- sign(gridy.temp) * gridy.temp^2
        }
        else if (gridtype1=="exp")
        {
            gridy.temp <- seq(exp(xmin), exp(xmax), length=gridsize)
            gridy <- log(gridy.temp)
        }
        gridtype.vec <- gridtype1
    }
    n <- length(y)

    est <- dnorm.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=w/n)
    fhat <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE)
  
    if (positive)
    {
        ## compute transformation KDE
        fhat$estimate <- fhat$estimate/(exp(gridy))
        fhat$x <- x
        fhat$eval.points <- exp(gridy) - adj.positive  
    }
    class(fhat) <- "kde"
  
    return(fhat)
}

kde.positive.1d <- function(x, h, adj.positive, binned=FALSE, xmin, xmax, compute.cont=TRUE, approx.cont=TRUE, ...)
{
    if (missing(adj.positive)) adj.positive <- abs(min(x)) 
    y <- log(x + adj.positive) 
    if (missing(h)) h <- hpi(y, binned=binned)

    d <- 1
    tol <- 3.7
    tol.h <-  tol*h
    if (missing(xmin)) xmin <- min(x) - tol.h
    if (missing(xmax)) xmax <- max(x) + tol.h
    xmin[xmin<0] <- 0
    ymin1 <- log(xmin + adj.positive)
    ymax1 <- log(xmax + adj.positive)

    fhaty <- kde(x=y, h=h, xmin=ymin1, xmax=ymax1, gridtype=c("exp"), binned=binned, compute.cont=compute.cont, approx.cont=approx.cont, ...) 
    fhaty$estimate[is.nan(fhaty$estimate)] <- 0

    fhatx <- fhaty
    fhatx$x <- x
    fhatx$eval.points <- exp(fhaty$eval.points) - adj.positive
    jacobian <- abs(exp(fhaty$eval.points))
    jacobian[jacobian<=0] <- min(fhatx$estimate[fhatx$estimate>0]) 
    fhatx$estimate <- fhaty$estimate/jacobian

    if (compute.cont)
        fhatx$cont <- contourLevels(fhatx, cont=1:99, approx=approx.cont)
      
    ## re-sample on regular grid    
    ep <- seq(fhatx$eval.points[1], tail(fhatx$eval.points,n=1), length=length(fhatx$eval.points))
    fhatx$estimate <- predict(fhatx, x=ep)
    fhatx$eval.points <- ep

    return(fhatx)
}

kde.unit.interval.1d <- function(x, h, binned=FALSE)
{
    d <- 1
    y <- qnorm(x)
    if (missing(h)) h <- hpi(y)
    xseq <- tail(head(seq(0,1, length=default.gridsize(d)+2),n=-1), n=-1)
    fhaty <- kde(x=y, h=h, binned=binned)
    fhaty$estimate <- predict(fhaty, x=qnorm(xseq))
    fhatx <- fhaty
    fhatx$eval.points <- xseq
    fhatx$estimate <- fhaty$estimate/dnorm(qnorm(xseq))

    ## apply loess smoothing for unsmooth binned estimates
    if (binned)
    {      
        fhatx.loess <- loess(fhatx$estimate ~ fhatx$eval.points)
        fhatx.smoothed <- fhatx
        fhatx.smoothed$eval.points <- xseq
        fhatx.smoothed$estimate <- predict(fhatx.loess, x=xseq)
        fhatx <- fhatx.smoothed
    }
    fhatx$x <- x

    return(fhatx)
}

###############################################################################
## Bivariate kernel density estimate using normal kernels, evaluated over grid
##
## Parameters
## x - data points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
###############################################################################

kde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, verbose=FALSE)
{
    ## initialise grid 
    n <- nrow(x)
    if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 

    suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

    if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)    
    fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
    if (verbose) pb <- txtProgressBar()
  
    for (i in 1:n)
    {
        ## compute evaluation points 
        eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
        eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
        eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
        eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
        eval.x.len <- length(eval.x)
        eval.pts <- expand.grid(eval.x, eval.y)
        fhat <- dmvnorm(eval.pts, x[i,], H)

        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        for (j in 1:length(eval.y))
        fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
        fhat.grid[eval.x.ind, eval.y.ind[j]] + 
          w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        if (verbose) setTxtProgressBar(pb, i/n)
    }
    if (verbose) close(pb)
    fhat.grid <- fhat.grid/n
    gridx1 <- list(gridx[[1]], gridx[[2]]) 
    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE)

    return(fhat.list)
}

######################################################################
## Bivariate KDE for data in positive quadrant
######################################################################

kde.positive.2d <- function(x, H, adj.positive, binned=FALSE, xmin, xmax, compute.cont=TRUE, approx.cont=TRUE, ...)
{
    if (missing(adj.positive)) adj.positive <- abs(apply(x, 2, min)) 
    y <- log(cbind(x[,1] + adj.positive[1],x[,2] + adj.positive[2])) 
    if (missing(H)) H <- Hpi(y, binned=binned)

    d <- ncol(x)
    tol <- 3.7
    tol.H <-  tol * diag(H)
    if (missing(xmin)) xmin <- apply(x, 2, min) - tol.H
    if (missing(xmax)) xmax <- apply(x, 2, max) + tol.H
    xmin[xmin<0] <- 0
    ymin1 <- log(pmax(xmin + adj.positive, apply(x, 2, min)))
    ymax1 <- log(xmax + adj.positive)

    fhaty <- kde(x=y, H=H, xmin=ymin1, xmax=ymax1, gridtype=c("exp", "exp"), binned=binned, compute.cont=compute.cont, approx.cont=approx.cont, ...) 
    fhaty$estimate[is.nan(fhaty$estimate)] <- 0

    fhatx <- fhaty
    fhatx$x <- x
    fhatx$eval.points[[1]] <- exp(fhaty$eval.points[[1]]) - adj.positive[1]
    fhatx$eval.points[[2]] <- exp(fhaty$eval.points[[2]]) - adj.positive[2]
    jacobian <- abs(exp(fhaty$eval.points[[1]]) %o% exp(fhaty$eval.points[[2]]))
    jacobian[jacobian<=0] <- min(fhatx$estimate[fhatx$estimate>0]) 
    fhatx$estimate <- fhaty$estimate/jacobian

    if (compute.cont)
        fhatx$cont <- contourLevels(fhatx, cont=1:99, approx=approx.cont)

    return(fhatx)
}

###############################################################################
## Trivariate kernel density estimate using normal kernels, evaluated over grid
##
## Parameters
## x - data points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
###############################################################################

kde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, verbose=FALSE)
{
    ## initialise grid 
    n <- nrow(x)

    if (is.null(gridx))
    gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
    suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

    if (is.null(grid.pts))
    grid.pts <- find.gridpts(gridx, suppx)    
    fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), length(gridx[[3]])))
    if (verbose) pb <- txtProgressBar()

    for (i in 1:n)
    {
        ## compute evaluation points
        eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
        eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
        eval.z <- gridx[[3]][grid.pts$xmin[i,3]:grid.pts$xmax[i,3]]
        eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
        eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
        eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
        eval.x.len <- length(eval.x)
        eval.pts <- expand.grid(eval.x, eval.y)

        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 

        for (k in 1:length(eval.z))
        {
            fhat <- w[i]*dmvnorm(cbind(eval.pts, eval.z[k]), x[i,], H)
            for (j in 1:length(eval.y))
                fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                    fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                    fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
        if (verbose) setTxtProgressBar(pb, i/n)
    }
    if (verbose) close(pb)
    fhat.grid <- fhat.grid/n
    gridx1 <- list(gridx[[1]], gridx[[2]], gridx[[3]]) 
    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE)

    return(fhat.list)
}

kde.grid.nd <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, verbose=FALSE)
{
    ## initialise grid 
    n <- nrow(x)
    if (is.null(gridx))
        gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
    
    gridx1 <- gridx
    gridx1$stepsize <- NULL
    gridx1$gridtype <- NULL
    eval.points <- do.call(expand.grid, gridx1)
    est <- kde.points(x=x, H=H, eval.points=eval.points, w=w, verbose=verbose)$estimate 
    fhat.grid <- array(est, dim=gridsize)
        
    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE)

    return(fhat.list)
}

###############################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
##
## Parameters
## x - data points
## H - bandwidth matrix
## eval.points - points where to evaluate density estimate
##
## Returns
## list with fields
## x - data points
## eval.points - points that KDE is evaluated at
## estimate - KDE evaluated at eval.points 
## H - bandwidth matrix 
###############################################################################

kde.points <- function(x, H, eval.points, w, verbose) 
{
    n <- nrow(x)
    d <- ncol(x)
    ne <- nrow(eval.points)
    Hs <- replicate(n, H, simplify=FALSE) 
    Hs <- do.call(rbind, Hs)
    fhat <- dmvnorm.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n, verbose=verbose)
    
    return(list(x=x, eval.points=eval.points, estimate=fhat, H=H, gridded=FALSE))
}

kde.points.1d <- function(x, h, eval.points, positive=FALSE, adj.positive, w) 
{
    n <- length(x)

    if (positive)
    {
        if (missing(adj.positive)) adj.positive <- abs(min(x))
        y <- log(x + adj.positive)  ## transform positive data x to real line
        eval.pointsy <- log(eval.points + adj.positive)
    }
    else
    {
        y <- x
        eval.pointsy <- eval.points
    }

    est <- dnorm.mixt(x=eval.pointsy, mus=y, sigmas=rep(h,n), props=w/n)
    if (positive)
        est <- est/(eval.points + adj.positive) 
    fhat <- list(x=x, eval.points=eval.points, estimate=est, h=h, H=h^2, gridded=FALSE)

    return(fhat)
}

#############################################################################
## S3 methods for KDE objects
#############################################################################

## predict method for KDE objects

predict.kde <- function(object, ..., x, zero.flag=TRUE)
{
    fhat <- grid.interp(x=x, gridx=object$eval.points, f=object$estimate)
    if (!zero.flag) warning("zero.flag=FALSE has been deprecated and no longer has any effect")
  
    return(fhat)
}

## plot method

plot.kde <- function(x, ...)
{ 
    fhat <- x
    if (is.vector(fhat$x)) plotkde.1d(fhat, ...)
    else
    {
        d <- ncol(fhat$x)

        if (d==2) 
        {
            opr <- options()$preferRaster; if (!is.null(opr)) if (!opr) options("preferRaster"=TRUE)
            plotret <- plotkde.2d(fhat, ...)
            if (!is.null(opr)) options("preferRaster"=opr)
            invisible(plotret)
        }
        else if (d==3)
        {
            plotkde.3d(fhat, ...)
            invisible()
        }
        else 
          stop ("Plot function only available for 1, 2 or 3-d data")
    }
}

plotkde.1d <- function(fhat, xlab, ylab="Density function", add=FALSE,
  drawpoints=FALSE, col=1, col.pt=4, col.cont=1, cont.lwd=1, jitter=FALSE, cont, abs.cont, approx.cont=TRUE, alpha=1, ...) 
{
    if (missing(xlab)) xlab <- fhat$names
    col <- transparency.col(col, alpha=alpha)

    if (add) lines(fhat$eval.points, fhat$estimate, xlab=xlab, ylab=ylab, col=col, ...)
    else plot(fhat$eval.points, fhat$estimate, type="l", xlab=xlab, ylab=ylab, col=col, ...) 

    ## compute contours
    if (!missing(cont) | !missing(abs.cont)) 
    {
        if (missing(abs.cont))
        {
            if (!is.null(fhat$cont))
            {
                cont.ind <- rep(FALSE, length(cont))
                for (j in 1:length(cont))
                {
                    ci <- which(cont[j]==(100-as.numeric(unlist(strsplit(names(fhat$cont),"%")))))
                    if (length(ci)==0) cont.ind[j] <- NA else cont.ind[j] = ci
                }
              
                if (any(is.na(cont.ind)))
                    hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
                else
                    hts <- fhat$cont[cont.ind]
            }
            else
                hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        }
        else
          hts <- abs.cont 
       
 
        if (is.null(fhat$deriv.order))
        {
            hts <- sort(hts, decreasing=TRUE)
            cont.ind <- 1-as.numeric(fhat$estimate>=hts[1])
            cont.ind[cont.ind==1] <- NA   
            lines(fhat$eval.points, cont.ind, col=col.cont, lwd=cont.lwd)
        }
        else
        {
           for (i in 1:length(hts))
           {     
               cont.ind <- 1-as.numeric(abs(fhat$estimate)>=abs(hts[i]))
               cont.ind[cont.ind==1] <- NA    
               lines(fhat$eval.points, cont.ind, col=col.cont, lwd=cont.lwd)      
           }
        }
    }
  
    if (drawpoints)
        if (jitter) rug(jitter(fhat$x), col=col.pt)
        else rug(fhat$x, col=col.pt)
}

plotkde.2d <- function(fhat, display="slice", cont=c(25,50,75), abs.cont, approx.cont=TRUE, xlab, ylab, zlab="Density function", cex=1, pch=1, labcex=1, add=FALSE, drawpoints=FALSE, drawlabels=TRUE, theta=-30, phi=40, d=4, col.pt=4, col, col.fun, alpha=1, lwd=1, border=1, thin=3, kdde.flag=FALSE, ticktype="detailed", ...) 
{
    disp <- match.arg(display, c("slice", "persp", "image", "filled.contour", "filled.contour2"))
    if (disp=="filled.contour2") disp <- "filled.contour"
    if (!is.list(fhat$eval.points)) stop("Need a grid of density estimates")

    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (missing(col.fun)) 
    {
        if (any(fhat$type=="kcurv")) col.fun <- function(n) { hcl.colors(n, palette="Oranges",rev=TRUE, alpha=alpha) }
        else col.fun <- function(n) { hcl.colors(n, palette="heat",rev=TRUE, alpha=alpha) }
    }
  
    ## perspective/wireframe plot
    if (disp=="persp")
    {
        hts <- seq(0, 1.1*max(fhat$estimate,na.rm=TRUE), length=100)
        if (missing(col)) col <- col.fun(length(hts)+1)
        if (length(col)<length(hts)) col <- rep(col, length=length(hts))
        col <- transparency.col(col, alpha=alpha)

        ## thinning indices
        plot.ind <- list(seq(1, length(fhat$eval.points[[1]]), by=thin), seq(1, length(fhat$eval.points[[2]]), by=thin))

        z <- fhat$estimate[plot.ind[[1]], plot.ind[[2]]]
        nrz <- nrow(z)
        ncz <- ncol(z)
        zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
        facetcol <- cut(zfacet, length(col)+1)

        plotret <- persp(fhat$eval.points[[1]][plot.ind[[1]]], fhat$eval.points[[2]][plot.ind[[2]]], z, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, col=col[facetcol], border=border, ticktype=ticktype, ...)
    }
    else if (disp=="slice") 
    {
        ## compute contours
        if (missing(abs.cont))
        {
            if (!is.null(fhat$cont))
            {
                cont.ind <- rep(FALSE, length(cont))
                for (j in 1:length(cont))
                {
                    ci <- which(cont[j]==(100-as.numeric(unlist(strsplit(names(fhat$cont),"%")))))
                    if (length(ci)==0) cont.ind[j] <- NA else cont.ind[j] = ci
                }
              
                if (any(is.na(cont.ind)))
                    hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
                else
                    hts <- fhat$cont[cont.ind]
            }
            else
                hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        }
        else
            hts <- abs.cont 

        hts <- sort(hts)
        if (missing(col)) col <- col.fun(length(hts))
        if (length(col)<length(hts)) col <- rep(col, times=length(hts))
        col <- transparency.col(col, alpha=alpha)

        ## draw contours
        j <- 0
        for (i in 1:length(hts)) 
        {
            if (missing(abs.cont)) { ni <- length(hts)-i+1; scale <- (100-cont[i])/hts[i]; scale2 <- cont[ni]/hts[ni] }
            else { ni <- i; scale <- 1; scale2 <-1 }

            if (hts[i]>0 | !is.null(fhat$deriv.order))
            {
                j <- j+1
                if (j==1) 
                    contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate*scale, level=hts[i]*scale, label=signif(hts[ni]*scale2), add=add, drawlabels=drawlabels, col=col[i], lwd=lwd, labcex=labcex, xlab=xlab, ylab=ylab, ...)
                else 
                    contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate*scale, level=hts[i]*scale, label=signif(hts[ni]*scale2), add=TRUE, drawlabels=drawlabels, col=col[i], lwd=lwd, labcex=labcex, ...)
            }
        }

        ## add points 
        if (drawpoints) points(fhat$x[,1], fhat$x[,2], col=col.pt, cex=cex, pch=pch)
    }
    ## image plot
    else if (disp=="image")
    {
        if (missing(col)) col <- col.fun(100)
        col <- transparency.col(col, alpha=alpha)
        image(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate, xlab=xlab, ylab=ylab, add=add, col=col, ...)

        ## add points 
        if (drawpoints) points(fhat$x[,1], fhat$x[,2], col=col.pt, cex=cex, pch=pch)
        box()
    }
    else if (disp=="filled.contour")
    {
        ## compute contours
        if (missing(abs.cont))
        {
            if (!is.null(fhat$cont))
            {
                cont.ind <- rep(FALSE, length(cont))
                for (j in 1:length(cont))
                {
                    ci <- which(cont[j]==(100-as.numeric(unlist(strsplit(names(fhat$cont),"%")))))
                    if (length(ci)==0) cont.ind[j] <- NA else cont.ind[j] = ci
                }
              
                if (any(is.na(cont.ind)))
                    hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
                else
                    hts <- fhat$cont[cont.ind]
            }
            else
                hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
        }
        else
            hts <- abs.cont 
        hts <- sort(hts)

        if (missing(col)) col <- c("transparent", col.fun(length(hts)))
        col <- transparency.col(col, alpha=alpha)
        clev <- c(min(c(fhat$estimate, hts)-0.01*max(abs(fhat$estimate))), hts, max(c(fhat$estimate, hts)) + 0.01*max(abs(fhat$estimate)))

        if (!add) plot(fhat$eval.points[[1]], fhat$eval.points[[2]], type="n", xlab=xlab, ylab=ylab, ...)
        .filled.contour(fhat$eval.points[[1]], fhat$eval.points[[2]], z=fhat$estimate, levels=clev, col=col)
      
        if (!missing(lwd))
        {
            for (i in 1:length(hts)) 
            {
                if (missing(abs.cont)) { ni <- length(hts)-i+1; scale <- (100-cont[i])/hts[i]; scale2 <- cont[ni]/hts[ni] }
                else { ni <- i; scale <- 1; scale2 <-1 }

                if (lwd >=1) contour(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$estimate*scale, level=hts[i]*scale, label=signif(hts[ni]*scale2,3), add=TRUE, drawlabels=drawlabels, col=1, lwd=lwd, labcex=labcex, ...)
            }
        }

        ## add points 
        if (drawpoints) points(fhat$x[,1], fhat$x[,2], col=col.pt, cex=cex, pch=pch)
        box()
    }

    if (disp=="persp") invisible(plotret)
    else invisible()
}
  
plotkde.3d <- function(fhat, display="plot3D", cont=c(25,50,75), abs.cont, approx.cont=TRUE, colors, col, col.fun, alphavec, size=3, cex=1, pch=1, theta=-30, phi=40, d=4, ticktype="detailed", bty="f", col.pt=4, add=FALSE, xlab, ylab, zlab, drawpoints=FALSE, alpha, box=TRUE, axes=TRUE, ...)
{
    ## compute contours
    if (missing(abs.cont))
    {
        if (!is.null(fhat$cont))
        {
            cont.ind <- rep(FALSE, length(fhat$cont))
            for (j in 1:length(cont))
                cont.ind[which(cont[j] == 100-as.numeric(unlist(strsplit(names(fhat$cont),"%"))))] <- TRUE
           
            if (all(!cont.ind))
                hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
            else
                hts <- fhat$cont[cont.ind]
         }
         else
            hts <- contourLevels(fhat, prob=(100-cont)/100, approx=approx.cont)
    }  
    else
        hts <- abs.cont
        
    nc <- length(hts)
    if (missing(col)) 
    { 
        if (missing(col.fun)) 
        {
            if (any(fhat$type=="kcurv")) col.fun <- function(n) { hcl.colors(n, palette="Oranges",rev=TRUE) }
            else col.fun <- function(n) { hcl.colors(n, palette="heat",rev=TRUE) }
        }
        col <- col.fun(n=length(hts))
    }    
    colors <- col
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (missing(zlab)) zlab <- fhat$names[3]
    if (missing(alphavec))
    {
        if (is.null(fhat$deriv.order)) alphavec <- seq(0.1,0.5,length=nc)
        else alphavec <- c(rev(seq(0.1,0.4,length=round(nc/2))), seq(0.1,0.4,length=round(nc/2)))
    }
    if (missing(alpha)) alpha <- 0.1 
    else if (!missing(alpha)) { alphavec <- rep(alpha,nc) }
            
    disp <- match.arg(display, c("plot3D", "rgl")) 
    if (disp %in% "plot3D")
    {
         
        for (i in 1:nc)
            if (hts[nc-i+1] < max(fhat$estimate))
                plot3D::isosurf3D(x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], colvar=fhat$estimate, level=hts[nc-i+1], add=add | (i>1), col=colors[i], alpha=alphavec[i], phi=phi, theta=theta, xlab=xlab, ylab=ylab, zlab=zlab, d=d, ticktype=ticktype, bty=bty, ...)
                
        if (drawpoints) plot3D::points3D(x=fhat$x[,1], y=fhat$x[,2], z=fhat$x[,3], cex=cex, col=col.pt, add=TRUE, pch=pch, d=d)    
    }
    else if (disp %in% "rgl")
    {   
        ## suggestions from Viktor Petukhov 08/03/2018
        if (!requireNamespace("rgl", quietly=TRUE)) stop("Install the rgl package as it is required.", call.=FALSE)
        if (!requireNamespace("misc3d", quietly=TRUE)) stop("Install the misc3d package as it is required.", call.=FALSE)
    
        fhat.eval.mean <- sapply(fhat$eval.points, mean)
        if (drawpoints)
            rgl::plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=size, col=col.pt, alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
        else
            rgl::plot3d(fhat$x[,1],fhat$x[,2],fhat$x[,3], size=0, col="transparent", alpha=0, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)

        for (i in 1:nc)
            if (hts[nc-i+1] < max(fhat$estimate))
                misc3d::contour3d(fhat$estimate, level=hts[nc-i+1], x=fhat$eval.points[[1]], y=fhat$eval.points[[2]], z=fhat$eval.points[[3]], add=TRUE, color=colors[i], alpha=alphavec[i], box=FALSE, axes=FALSE, ...)
    
        if (axes) rgl::axes3d(c("x","y","z"))
        if (box) rgl::box3d()
    }
}

## contourLevels method
## create S3 generic 

contourLevels <- function(x, ...) { UseMethod("contourLevels") }   

contourLevels.kde <- function(x, prob, cont, nlevels=5, approx=TRUE, ...)
{ 
    fhat <- x
    if (is.vector(fhat$x))
    {
        d <- 1; n <- length(fhat$x)
    }
    else
    {
        d <- ncol(fhat$x); n <-nrow(fhat$x)
        if (!is.matrix(fhat$x)) fhat$x <- as.matrix(fhat$x)
    }

    if (is.null(x$w)) w <- rep(1, n)
    else w <- x$w

    if (is.null(fhat$gridded))
    {
        if (d==1) fhat$gridded <- fhat$binned
        else fhat$gridded <- is.list(fhat$eval.points)
    }
  
    if (missing(prob) & missing(cont))
        hts <- pretty(fhat$estimate, n=nlevels) 
    else
    {
        if (approx & fhat$gridded)
            dobs <- predict(fhat, x=fhat$x)
        else
            dobs <- kde(x=fhat$x, H=fhat$H, eval.points=fhat$x, w=w)$estimate 
        
        if (!missing(prob) & missing(cont))
            hts <- quantile(dobs, prob=prob)
          
        if (missing(prob) & !missing(cont))
            hts <- quantile(dobs, prob=(100-cont)/100)
    }
      
    return(hts)
}

###############################################################################
## Riemann sums to compute approximate Lebesgue measure of contour set
###############################################################################

contourSizes <- function(x, abs.cont, cont=c(25,50,75), approx=TRUE)
{
    if (missing(abs.cont))
        abs.cont <- contourLevels(x, cont=cont, approx=approx)
    num.int <- rep(0, length(abs.cont))
    if (!is.null(names(abs.cont))) names(num.int) <- names(abs.cont) 
    
    if (!is.list(x$eval.points))
        delta.int <- head(diff(x$eval.points), n=1)
    else
        delta.int <- prod(sapply(lapply(x$eval.points, diff), head, n=1)) 
    
    for (j in 1:length(abs.cont)) 
        num.int[j] <- sum(x$estimate>abs.cont[j])
    
    return(num.int*delta.int)
}

###############################################################################
## Riemann sums to compute approximate probability of contour set
###############################################################################

contourProbs <- function(x, abs.cont, cont=c(25,50,75), approx=TRUE)
{
    if (missing(abs.cont)) 
        abs.cont <- contourLevels(x, cont=cont, approx=approx)
    num.int <- rep(0, length(abs.cont))
    if (!is.null(names(abs.cont))) 
        names(num.int) <- names(abs.cont)
    if (!is.list(x$eval.points)) 
    {
        delta.int <- head(diff(x$eval.points), n=1)
        eval.points.midpoint <- (head(x$eval.points,n=-1)+tail(x$eval.points,n=-1))/2
    }
    else
    {
        delta.int <- prod(sapply(lapply(x$eval.points, diff), head, n=1)) 
        eval.points.midpoint <- expand.grid(lapply(x$eval.points, function(y) { (head(y,n=-1)+tail(y,n=-1))/2 }))
    }
    
    x.evmp <- predict(x, x=eval.points.midpoint)
    for (i in 1:length(num.int)) num.int[i] <- sum(x.evmp*(x.evmp>=abs.cont[i])*delta.int)
    
    return(num.int)
}
 

###############################################################################
## Generate grid over a set of points
##
## Parameters
## x - data points
## H - bandwidth matrix
## tol - tolerance = extra coverage exceeding the range of x   
## gridsize - number of points for each direction
##
## Returns
## gridx - list of intervals, one for each co-ord direction so that
##         gridx[[1]] x gridx[[2]] x ... x gridx[[d]] is the grid
## stepsize - vector of step sizes 
###############################################################################

make.grid.ks <- function(x, H, tol, gridsize, xmin, xmax, gridtype)
{
    d <- ncol(x)
    tol.H <-  tol * diag(H)
    if (missing(xmin)) xmin <- apply(x, 2, min) - tol.H
    if (missing(xmax)) xmax <- apply(x, 2, max) + tol.H

    stepsize <- rep(0, d)
    gridx <- numeric(0)
    if (length(gridsize)==1)  gridsize <- rep(gridsize, d)
    if (missing(gridtype)) gridtype <- rep("linear", d)
    gridtype.vec <- rep("", d)

    for (i in 1:d)
    {
        gridtype1 <- match.arg(gridtype[i], c("linear", "sqrt", "quantile", "exp")) 
        if (gridtype1=="linear")
        {  
            gridx <- c(gridx, list(seq(xmin[i], xmax[i], length=gridsize[i])))
            stepsize[i] <- abs(gridx[[i]][1] - gridx[[i]][2])
        }
        else if (gridtype1=="sqrt")
        {
            gridx.temp <- seq(sign(xmin[i])*sqrt(abs(xmin[i])), sign(xmax[i])*sqrt(abs(xmax[i])), length=gridsize[i])
            gridx <- c(gridx, list(sign(gridx.temp) * gridx.temp^2))
            stepsize[i] <- NA
            }
        else if (gridtype1=="quantile")
        {
            gridx.temp <- qnorm(seq(1e-2, 1-1e-2, length=gridsize[i]))
            gridx <- c(gridx, list(xmin[i] + (xmax[i]-xmin[i])*(gridx.temp-min(gridx.temp))/(max(gridx.temp)-min(gridx.temp))))
            stepsize[i] <- NA
        }
        else if (gridtype1=="exp")
        {
            gridx.temp <- seq(exp(xmin[i]), exp(xmax[i]), length=gridsize[i])
            gridx <- c(gridx, list(log(gridx.temp)))
            stepsize[i] <- NA
        }
        gridtype.vec[i] <- gridtype1
    }
    gridx <- c(gridx, list(stepsize = stepsize, gridtype=gridtype.vec))

    return(gridx)
}  

###############################################################################
## Generate kernel (rectangular) support at data point
## 
## Parameters
## x - data points
## H - bandwidth matrix
## tol - tolerance = extra coverage exceeding the range of x 
##
## Returns
## list of min and max points of support (here we parameterise rectangles
## by their min = lower left co-ord and max = upper right coord)
###############################################################################

make.supp <- function(x, H, tol)
{
    n <- nrow(x)
    d <- ncol(x)
    tol.H <- tol * diag(H)
    xmin <- matrix(0, nrow=n, ncol=d)
    xmax <- matrix(0, nrow=n, ncol=d)

    for (i in 1:n)
    {
        xmin[i,] <- x[i,] - tol.H
        xmax[i,] <- x[i,] + tol.H 
    }
           
    return(list(xmin=xmin, xmax=xmax))
}

###############################################################################
## Find the grid points contained in kernel support rectangles 
##
## Parameters
## gridx - grid (list of subdivided intervals)
## rectx - rectangles (list of min and max points) 
##
## Returns
## list of min and max points of the grid for each rectangle 
###############################################################################

find.gridpts <- function(gridx, suppx)
{
    xmax <- suppx$xmax
    xmin <- suppx$xmin
    d <- ncol(xmax)
    n <- nrow(xmax)
    gridpts.min <- matrix(0, ncol=d, nrow=n)
    gridpts.max <- gridpts.min
  
    for (i in 1:n)
        for (j in 1:d)    
        {
            ## find index of last element of gridx smaller than min support  
            tsum <- sum(xmin[i,j] >= gridx[[j]])
            if (tsum==0)
                gridpts.min[i,j] <- 1
            else
                gridpts.min[i,j] <- tsum

            ## find index of first element gridx greater than max support 
            gridpts.max[i,j] <- sum(xmax[i,j] >= gridx[[j]])
        }   
        
    return(list(xmin=gridpts.min, xmax=gridpts.max))
} 

##############################################################################
## Interpolate the values of f defined on gridx at new values x 
##############################################################################

grid.interp <- function(x, gridx, f)
{
    if (!is.list(gridx))
    {
        ## uniform grid
        if (isTRUE(all.equal(diff(gridx), rep(diff(gridx)[1], length(gridx)-1))))
            fx <- grid.interp.1d(x=as.vector(x), gridx=gridx, f=f)
        ## non-uniform grid
        else
            fx <- varying.grid.interp.1d(x=as.vector(x), gridx=gridx, f=f)
    
    }
    else
    {  
        if (is.vector(x)) x <- as.matrix(t(x))
        d <- ncol(x)
        n <- nrow(x)

        if (d<2) stop("x should be a vector") 
        gridx.diff <- lapply(lapply(gridx,diff), getElement, 1)
        for (i in 1:length(gridx.diff)) gridx.diff[[i]] <- rep(gridx.diff[[i]], sapply(gridx, length)[i]-1)
        uniform.grid.flag <- isTRUE(all.equal(lapply(gridx,diff), gridx.diff))

        ## uniform grid
        if (d==2 & uniform.grid.flag) fx <- grid.interp.2d(x=x, gridx=gridx, f=f)
        else if (d==3 & uniform.grid.flag) fx <- grid.interp.3d(x=x, gridx=gridx, f=f)
        else ## d >=4 or non-uniform grid
        {    
            gridsize <- sapply(gridx,length)
            gind <- matrix(0, nrow=n, ncol=d)
            
            for (i in 1:n)
                for (j in 1:d)
                {
                    tsum <- sum(x[i,j] >= gridx[[j]])
                    if (tsum==0) gind[i,j] <- 1
                    else gind[i,j] <- tsum
                }
            for (j in 1:d) gind[gind[,j]>=gridsize[j],j] <- gridsize[j]-1

            bperm <- list()
            for (j in 1:d) bperm[[j]] <- elem(1,2)
            binary.perm <- as.matrix(expand.grid(bperm))
            colnames(binary.perm) <- NULL
            
            gind.list <- list()
            fx <- rep(0, length=n)
            for (i in 1:n)
            {
                gind.list[[i]] <- matrix(gind[i,], nrow=2^d, ncol=d, byrow=TRUE) + binary.perm
                w <- matrix(0, nrow=2^d, ncol=d)
                gridw <- matrix(0, nrow=2^d, ncol=d)
                for (j in 1:d)
                {
                    gind.list[[i]][,j][gind.list[[i]][,j]>=gridsize[j]] <- gridsize[j]
                    gridw[,j] <- gridx[[j]][gind.list[[i]][,j]]
                }
                w <- abs(matrix(as.numeric(x[i,]), nrow=2^d, ncol=d, byrow=TRUE) - gridw)
                w <- apply(w, 1, prod)
                w <- w/sum(w)
                fx[i] <- sum(w*f[gind.list[[i]][2^d:1,]])
            }
        }
    }
  
    return(fx)
}

grid.interp.1d <- function(x, gridx, f)
{
   n <- length(x)
   gpoints1 <- gridx
   M1 <- length(gpoints1)
   a1 <- gpoints1[1]
   b1 <- gpoints1[M1]
   
   out <- .C(C_interp1d, x1=as.double(x), n=as.integer(n), a1=as.double(a1), b1=as.double(b1), M1=as.integer(M1), fun=as.double(as.vector(f)), est=double(n))
   return(out$est)
}

grid.interp.2d <- function(x, gridx, f)
{
   n <- nrow(x)
   gpoints1 <- gridx[[1]]
   gpoints2 <- gridx[[2]]
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]

   out <- .C(C_interp2d, x1=as.double(x[,1]), x2=as.double(x[,2]), n=as.integer(n), a1=as.double(a1), a2=as.double(a2), b1=as.double(b1), b2=as.double(b2), M1=as.integer(M1), M2=as.integer(M2), fun=as.double(as.vector(f)), est=double(n))
   
   return(out$est)
}

grid.interp.3d <- function(x, gridx, f)
{
   n <- nrow(x)
   gpoints1 <- gridx[[1]]
   gpoints2 <- gridx[[2]]
   gpoints3 <- gridx[[3]]
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   
   out <- .C(C_interp3d, x1=as.double(x[,1]), x2=as.double(x[,2]), x3=as.double(x[,3]), n=as.integer(n), a1=as.double(a1), a2=as.double(a2), a3=as.double(a3), b1=as.double(b1), b2=as.double(b2), b3=as.double(b3), M1=as.integer(M1), M2=as.integer(M2), M3=as.integer(M3), fun=as.double(as.vector(f)), est=double(n))
   
   return(out$est)
}

## Linear intepolation based on kernel estimation grid
## alias for predict.kde

kde.approx <- function(fhat, x)
{
    return(grid.interp(x=x, gridx=fhat$eval.points, f=fhat$estimate))
}

##############################################################################
## Find the nearest grid points surrounding point x for non-uniform grids 
##############################################################################

varying.grid.interp <- function(x, gridx, f)
{
    if (!is.list(gridx))
        return(varying.grid.interp.1d(x=x, gridx=gridx, f=f))
    else
    {  
        if (is.vector(x)) x <- as.matrix(t(x))
        d <- ncol(x)
        n <- nrow(x)
        gridsize <- sapply(gridx,length)
        gind <- matrix(0, nrow=n, ncol=d)

        for (i in 1:n)
            for (j in 1:d)
            {
                tsum <- sum(x[i,j] >= gridx[[j]])
                if (tsum==0) gind[i,j] <- 1
                else gind[i,j] <- tsum
            }
    }

    bperm <- list()
    for (j in 1:d) bperm[[j]] <- elem(1,2)
    binary.perm <- as.matrix(expand.grid(bperm))
    colnames(binary.perm) <- NULL

    gind.list <- list()
    fx <- rep(0, length=n)
    for (i in 1:n)
    {
        gind.list[[i]] <- matrix(gind[i,], nrow=2^d, ncol=d, byrow=TRUE) + binary.perm
        w <- matrix(0, nrow=2^d, ncol=d)
        gridw <- matrix(0, nrow=2^d, ncol=d)
        for (j in 1:d)
        {
            gind.list[[i]][,j][gind.list[[i]][,j]>=gridsize[j]] <- gridsize[j]
            gridw[,j] <- gridx[[j]][gind.list[[i]][,j]]
        }
        w <- 1/apply((matrix(as.numeric(x[i,]), nrow=2^d, ncol=d, byrow=TRUE) - gridw)^2, 1, sum)
        w[w>1e5] <- 1e5
        w <- w/sum(w)
        fx[i] <- sum(w*f[gind.list[[i]]])
    }

    return(fx)
}

varying.grid.interp.1d <- function(x, gridx, f)
{
    n <- length(x)
    gind <- rep(0, length=n)

    for (i in 1:length(x))
    {
        tsum <- sum(x[i] >= gridx)
        if (tsum==0) gind[i] <- 1
        else gind[i] <- tsum
    }
 
    gind2 <- gind+1
    gind2[gind2>length(gridx)] <- length(gridx)
    gind2[x<=gridx[1]] <- gind[x<=gridx[1]]
    gind <- cbind(gind, gind2)
    colnames(gind) <- NULL

    fx <- rep(0, n)
    for (i in 1:n)
    {
        w <- 1/(x[i] - gridx[gind[i,]])^2
        w[w>1e5] <- 1e5
        w <- w/sum(w)   
        fx[i] <- sum(w*f[gind[i,]])
    }

  return(fx)
}
