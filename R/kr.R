###############################################################################
## Multivariate kernel regression
###############################################################################

kr <- function(x, y, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE, xsupp.truncate=TRUE, reg.order=1)
{
    ## default values 
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize

    if (is.data.frame(x)) x <- as.matrix(x)
    p <- reg.order
    if (d==1 & missing(h)) h <- hpi(x=x, binned=default.bflag(d=d, n=n))  
    if (d>1 & missing(H) & d>1) H <- Hpi(x=x, binned=default.bflag(d=d, n=n))
 
    ## compute binned estimator
    if (binned)
    {
        stop("Binning not yet available.")
    }
    else
    {
        ## compute exact (non-binned) estimator
       
        ## 1-dimensional    
        if (d==1)
        {
            ##if (missing(eval.points))
            ##{
            ##    if (unit.interval)
            ##    {
            ##        mhat <- kr.unit.interval.1d(x=x, y=y, h=h, binned=FALSE, reg.order=p)
            ##    }
            ##    else
            ##        mhat <- kr.grid.1d(x=x, y=y, h=h, gridsize=gridsize, supp=supp, positive=positive, xmin=xmin, xmax=xmax, adj.positive=adj.positive, gridtype=gridtype, w=w, reg.order=p)
           ## }
           ## else
           ##     mhat <- kr.points.1d(x=x, y=y, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive, w=w, reg.order=p)
        }
        ## multi-dimensional
        else
        {  
            if (missing(eval.points))
            {
                if (d==2)
                    mhat <- kr.grid.2d(x=x, y=y, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose, xsupp.truncate=xsupp.truncate, reg.order=p)
                ##else if (d == 3)
                ##    mhat <- kr.grid.3d(x=x, y=y, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose, reg.order=p) 
                else 
                    stop("Need to specify eval.points for more than 3 dimensions")
            }
            ##else
            ##    mhat <- kr.points(x=x, y=y, H=H, eval.points=eval.points, w=w, reg.order=p)     
        }
    }

    mhat$binned <- binned
    mhat$names <- parse.name(x)  ## add variable names
    mhat$w <- w
    class(mhat) <- "kr"
    
    ## compute prob contour levels
    if (compute.cont & missing(eval.points))
    {
        mhat.temp <- mhat
        mhat.temp$estimate[is.na(mhat.temp$estimate)] <- -1e10
        mhat$cont <- contourLevels(mhat.temp, cont=1:99, approx=approx.cont)
    }
    
    return(mhat)
}


kr.1d <- function(x, y, h, eval.points, reg.order=1)
{
    n <- length(x)
    nsumy <- sum(y)/n
    wy <- y/nsumy
    if (missing(eval.points))
        s0 <- kde(x=x, h=h)
    else
        s0 <- kde(x=x, h=h, eval.points=eval.points)
    
    fhat.wy <- kde(x=x, h=h, w=wy, eval.points=s0$eval.points)
    fhat.wy$estimate <- fhat.wy$estimate * nsumy
    mhat <- fhat.wy
    if (reg.order==0)
    {  
        mhat$estimate <- mhat$estimate/s0$estimate
    }
    
    if (reg.order==1)
    {
        fhat1.wy <- kdde(x=x, h=h, w=wy, eval.points=s0$eval.points, deriv.order=1)
        fhat1.wy$estimate <- fhat1.wy$estimate * nsumy*h
        s1 <- -kdde(x=x, h=h, deriv.order=1, eval.points=s0$eval.points)$estimate*h
        s2 <- kdde(x=x, h=h, deriv.order=2, eval.points=s0$eval.points)$estimate*h^2 - s0$estimate
        
        mhat$estimate <- s2/(s2*s0$estimate-s1^2) * fhat.wy$estimate + s1/(s2*s0$estimate-s1^2) * fhat1.wy$estimate
    }
    
    return(mhat)
}

kr.grid.2d <- function(x, y, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, verbose=FALSE, xsupp.truncate=FALSE, tol.zero=1e-3, reg.order=1)
{
   p <- reg.order

   ## initialise grid 
   n <- nrow(x)
   if (is.null(gridx))
       gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 

   mhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
   if (verbose) pb <- txtProgressBar()

   eval.pts <- expand.grid(gridx[[1]], gridx[[2]])
   eval.len <- nrow(eval.pts)

   if (xsupp.truncate)
   {
       fhat <- kde(x=x, binned=TRUE)
       fhat.cont <- contourLevels(fhat, cont=99, approx=TRUE)
       fhat.eval <- predict(fhat, x=eval.pts)
   }
   else
   {
       fhat.cont <- 0
       fhat.eval <- rep(1, eval.len)
   }
   
   beta.hat <- rep(NA, length=eval.len)
   
   for (i in 1:eval.len)
   {
       if (fhat.eval[i]>fhat.cont)
       {
           eval.ptsi <- unlist(eval.pts[i,])
           if (p==0) design.x <- rep(1, n)
           else if (p==1) design.x <- cbind(Inter=rep(1,n), sweep(x, 2, eval.ptsi))
           W.x <- dmvnorm(x=x, mean=eval.ptsi, sigma=H)
           A <- t(design.x * W.x)
           b <- A %*% y
           A <- A %*% design.x
           detA <- det(A)
           if (!(is.infinite(detA) | is.na(detA)))
               if (detA > tol.zero) beta.hat[i] <- solve(A,b)[1]
       }
       if (verbose) setTxtProgressBar(pb, i/eval.len)
   }
   
   mhat.grid <- array(beta.hat, dim=gridsize)
   if (verbose) close(pb)
   gridx1 <- list(gridx[[1]], gridx[[2]]) 
   
   mhat.list <- list(x=x, y=y, eval.points=gridx1, estimate=mhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, reg.order=p)
   
   return(mhat.list)
}


plot.kr <- function(x, display="rgl", col, col.fun, xlab, ylab, ...)
{
    if (display!="rgl") plot.kdde(x=x, col=col, col.fun=col.fun, display=display, xlab=xlab, ylab=ylab, ...)
    else
    {
        if (!requireNamespace("rgl", quietly=TRUE)) stop("Install the rgl package as it is required.", call.=FALSE)
        
        if (missing(col.fun)) col.fun <- terrain.colors
        if (missing(col)) col <- terrain.colors(5)[2]
        col.table <- col.fun(round(diff(range(x$estimate, na.rm=TRUE)) + 1,0))
        col <- col.table[x$estimate - min(x$estimate, na.rm=TRUE) + 1]
        if (missing(xlab)) xlab <- x$names[1]
        if (missing(ylab)) ylab <- x$names[2]
        rgl::persp3d(x$eval.points[[1]], x$eval.points[[2]], x$estimate, col=col, xlab=xlab, ylab=ylab, ...)
    }
}
   
contourLevels.kr <- function(x, ...)
{
    mhat <- x
    mhat$estimate[is.na(mhat$estimate)] <- 0
    mhat$deriv.order <- 1
    mhat$deriv.ind <- matrix(1, ncol=1, nrow=1)
    mhat$estimate <- list(mhat$estimate) 
    return(contourLevels.kdde(x=mhat,...))
}     

predict.kr <- function(object, ..., x)
{
    mhat <- object
    mhat$estimate[is.na(mhat$estimate)] <- 0
    return(predict.kde(object=mhat, ..., x=x))
}
