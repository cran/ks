######################################################################
## Balloon variable KDE
######################################################################

kde.balloon <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE)
{
    ## default values
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize

    ## clip data to xmin, xmax grid limits for binned estimation
    grid.clip <- binned    
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }

    ## default bandwidths
    if (d==1 & missing(h)) h <- hns(x=x, deriv.order=2)
    if (d>1 & missing(H)) H <- Hns(x=x, deriv.order=2) 
    
    if (d==2) fhat <- kde.balloon.2d(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont, verbose=verbose)
    else stop("kde.balloon only implemented for d=2")

    if (compute.cont)
        fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)

    return(fhat)
}

######################################################################
## Bivariate balloon variable KDE
######################################################################

kde.balloon.2d <- function(x, H, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE)
{
    d <- ncol(x)
    n <- nrow(x)
    fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont)
    fhat.ep <- expand.grid(fhat$eval.points)
   
    fhat.pilot <- fhat
    fhat2.pilot <- kdde(x=x, deriv.order=2, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w)

    h.pi <- (d*predict(fhat.pilot, x=fhat.ep)/drop((4*pi)^(d/2)*predict(fhat2.pilot, x=fhat.ep) %*% vec(diag(d)))^2)^(1/(d+4))*n^(-1/(d+4))
    h.pi <- array(h.pi, dim=dim(fhat$estimate))
    gs <- dim(fhat$estimate)
    
    if (verbose) { pb <- txtProgressBar(max=prod(gs)); k <- 0 }

    fhat$estimate <- array(0, dim=dim(fhat$estimate))
    for (i in 1:gs[2])
        for (j in 1:gs[1])
        {
            if (!is.na(h.pi[i,j]) & !is.infinite(h.pi[i,j])) if (h.pi[i,j]>0) fhat$estimate[i,j] <- kde(x=x, w=w, H=h.pi[i,j]^2*diag(d), eval.points=fhat.ep[i+(j-1)*gs[1],])$estimate
            if (verbose) { k <- k+1; setTxtProgressBar(pb,k) }
        }
    if (verbose) close(pb)

    ## re-scale density estimate to integral 1
    delta.int <- prod(sapply(fhat$eval.points, diff)[1,]) 
    riemann.sum <- sum(fhat$estimate*delta.int)
    fhat$estimate <- fhat$estimate/riemann.sum
    
    fhat$names <- parse.name(x)  ## add variable names
    if (compute.cont)
        fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)

    fhat$H <- h.pi^2
    
    return(fhat)
}

######################################################################
## Sample point variable KDE
######################################################################

kde.sp <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE)
{
    ## default values
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize

    ## clip data to xmin, xmax grid limits for binned estimation
    grid.clip <- binned    
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }

    ## default bandwidths
    if (d==1 & missing(h)) h <- hns(x=x, deriv.order=4)
    if (d>1 & missing(H)) H <- Hns(x=x, deriv.order=4) 

    if (d==2) fhat <- kde.sp.2d(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont, verbose=verbose, pre=TRUE)
    else stop("kde.sp only implemented for d=2")

    fhat$names <- parse.name(x)  ## add variable names
    if (compute.cont)
        fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)

    return(fhat)
}

######################################################################
## Bivariate sample point variable KDE
######################################################################

kde.sp.2d <- function(x, H, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE, pre=TRUE)
{
    d <- 2; n <- nrow(x)
    if (pre) 
    { 
        x.orig <- x
        S12 <- diag(apply(x.orig, 2, sd))  
        Sinv12 <- matrix.pow(S12,-1)
        x <- pre.scale(x)
        if (!missing(xmin)) xmin <- xmin %*% Sinv12
        if (!missing(xmax)) xmax <- xmax %*% Sinv12
        H <- Hns(x=x, deriv.order=4)
    }
        
    fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont)
    fhat.pilot <- fhat
    fhat1.pilot <- kdde(x=x, deriv.order=1, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w)
    fhat2.pilot <- kdde(x=x, deriv.order=2, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w)
    fhat3.pilot <- kdde(x=x, deriv.order=3, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w)
    fhat4.pilot <- kdde(x=x, deriv.order=4, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w)

    fhat.pilot <- predict(fhat.pilot, x=x)
    fhat1.pilot <- predict(fhat1.pilot, x=x)
    fhat2.pilot <- predict(fhat2.pilot, x=x)
    fhat3.pilot <- predict(fhat3.pilot, x=x)
    fhat4.pilot <- predict(fhat4.pilot, x=x)

    lambda1 <- 1/fhat.pilot^5*rowKpow(fhat1.pilot, r=4)
    lambda2 <- 1/fhat.pilot^4*rowKpow(fhat2.pilot, fhat1.pilot, r=1, s=2)
    lambda3 <- 1/fhat.pilot^3*rowKpow(fhat2.pilot, r=2)
    lambda4 <- 1/fhat.pilot^3*rowKpow(fhat3.pilot, fhat1.pilot, r=1, s=1)
    lambda5 <- 1/fhat.pilot^2*fhat4.pilot
    lambda <- drop((24*lambda1 - 36*lambda2 + 6*lambda3 + 8*lambda4 - lambda5) %*% Sdr(d,r=4) %*% (vec(diag(d)) %x% vec(diag(d))))

    RK <- (4*pi)^(-d/2)
    h.Ab  <- (8*d*RK*fhat.pilot^(1+d/2)/lambda^2)^(1/(d+8))*n^(-1/(d+8))
    
    fhat$estimate <- array(0, dim=dim(fhat$estimate))
    xmin <- sapply(fhat$eval.points, min)
    xmax <- sapply(fhat$eval.points, max)
    
    if (verbose) { pb <- txtProgressBar(max=n) }
    for (i in 1:n)
    {
        if (verbose) setTxtProgressBar(pb, i)
        HAb <- h.Ab[i]^2*diag(d)
        fhat$estimate <- fhat$estimate + kde(x=matrix(x[i,], nrow=1), H=HAb, xmin=xmin, xmax=xmax, binned=binned, gridsize=dim(fhat$estimate), bgridsize=dim(fhat$estimate))$estimate
    }
    if (verbose) close(pb)
    
    fhat$estimate <- fhat$estimate/n
    if (pre)
    {
        ep <- cbind(fhat$eval.points[[1]], fhat$eval.points[[2]]) %*% S12
        fhat$eval.points[[1]] <- ep[,1]
        fhat$eval.points[[2]] <- ep[,2]
        fhat$estimate <- fhat$estimate/det(S12)
        fhat$x <- x.orig
    }
    fhat$cont <- contourLevels(fhat, cont=1:99)
    fhat$H <- h.Ab^2
    
    return(fhat)
}

