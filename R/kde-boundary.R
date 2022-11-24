######################################################################
## Boundary KDE
######################################################################

kde.boundary <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, boundary.supp, boundary.kernel="beta", verbose=FALSE)
{
    bk <- match.arg(boundary.kernel, c("beta", "linear")) 
    if (bk=="beta")
    {
        if (missing(boundary.supp)) boundary.supp <- 10
        fhat <- kde.beta.boundary(x=x, H=H, h=h, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont, boundary.supp=boundary.supp, verbose=verbose) 
      
    }
    else if (bk=="linear")
    {
        if (missing(boundary.supp)) boundary.supp <- 2
        fhat <- kde.linear.boundary(x=x, H=H, h=h, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, compute.cont=compute.cont, approx.cont=approx.cont, boundary.supp=boundary.supp, verbose=verbose)
    }

    return(fhat)
}

######################################################################
## Linear boundary KDE
######################################################################

kde.linear.boundary <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, boundary.supp=2, verbose=FALSE)
{
    ## default values 
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    if (missing(binned)) binned <- ksd$binned
    if (missing(gridsize)) gridsize <- ksd$gridsize
    if (missing(bgridsize)) bgridsize <- gridsize
  
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
    if (d==1 & missing(h)) 
        h <- hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n))
    if (d>1 & missing(H))
        H <- Hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n))
    
    ## compute exact (non-binned) estimator
    ## 1-dimensional    
    if (d==1)
    {
        stop("Not yet implemented.")  
    }
    ## multi-dimensional
    else
    {  
        if (is.data.frame(x)) x <- as.matrix(x)
        
        if (missing(eval.points))
        {
            if (d==2)
                fhat <- kde.LB.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose)
            else 
                stop("Not yet implemented.")
        }
        else
            stop("Not yet implemented.")
    }
    
    fhat$binned <- binned
    fhat$names <- parse.name(x)  ## add variable names
    fhat$w <- w
    class(fhat) <- "kde"
    
    ## compute prob contour levels
    if (compute.cont & missing(eval.points))
        fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)
    
    return(fhat)
}

######################################################################
## Bivariate linear boundary KDE
######################################################################

kde.LB.grid.2d <- function(x, H, gridsize, bgridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, boundary.supp=10, binned=FALSE, verbose=FALSE)
{
   n <- nrow(x)
   d <- ncol(x)
   
   if (missing(xmin)) xmin <- apply(x, 2, min) 
   if (missing(xmax)) xmax <- apply(x, 2, max) 
   if (missing(gridtype)) gridtype <- rep("linear", d)
   
   h <- sqrt(diag(H))

   ## initialise grid 
   if (is.null(gridx))
       gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
   suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
   
   if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)
   fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
   
   ## indicator for closeness to boundary 
   bound.ind <- boundary.ind(x=x, h=h, boundary.supp=boundary.supp)
   n1 <- sum(!bound.ind)

   if (verbose) pb <- txtProgressBar()
   if (binned)
   {
      ## interior points - use normal kernel
      fhat.grid <- n1*kde(x=x[!bound.ind,], H=H, xmin=xmin, xmax=xmax, binned=binned, bgridsize=bgridsize)$estimate
       
        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)
            if (bound.ind[i])
            {
                ## compute evaluation points 
                eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
                eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
                eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
                eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
                eval.x.len <- length(eval.x)
                 
                ## use linear boundary kernel for boundary points 
                fhat <- dmvnorm.LB(x=expand.grid(eval.x, eval.y), mu=x[i,], Sigma=H)
             
                ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
                for (j in 1:length(eval.y))
                    fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
                        fhat.grid[eval.x.ind, eval.y.ind[j]] + 
                        w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
            }
        }
    }
    else
    {  
        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)

            ## compute evaluation points 
            eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
            eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
            eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
            eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
            eval.x.len <- length(eval.x)

            ## interior points - use normal kernel
            if (!bound.ind[i])
            {
                eval.pts <- expand.grid(list(eval.x, eval.y))
                fhat <- dmvnorm.mixt(x=eval.pts, mus=x[i,], Sigmas=H, props=1)
            }
            else
            {
                ## use linear boundary kernel for boundary points 
                 fhat <- dmvnorm.LB(x=expand.grid(eval.x, eval.y), mu=x[i,], Sigma=H)               
            }
             
            ## place vector of density estimate values fhat onto grid fhat.grid
            for (j in 1:length(eval.y))
                fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
                    fhat.grid[eval.x.ind, eval.y.ind[j]] + 
                    w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
    }
    fhat.grid <- fhat.grid/n
    gridx1 <- list(gridx[[1]], gridx[[2]])
    fhat.grid <- fhat.grid/sum(fhat.grid*prod(sapply(gridx1,diff)[1,]))

    if (verbose) close(pb)
    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, boundary=bound.ind)
   
    return(fhat.list)
}

## bivariate linear boundary normal kernel 

dmvnorm.LB.kernel.2d <- function(x, H, xmin=c(0,0), xmax=c(1,1), ...)
{
    x1 <- seq(xmin[1], xmax[1], length=151)
    x2 <- seq(xmin[2], xmax[2], length=151)
    eval.points <- list(x1, x2)

    fhat <- list()
    fhat$eval.points <- eval.points
    fhat$estimate <- array(dmvnorm.LB(x=expand.grid(eval.points), mu=x, Sigma=H), dim=c(length(x1), length(x2)))

    x <- rmvnorm.mixt(n=1000, mus=x, Sigmas=H, props=1) 
    fhat$x <- x
    fhat$H <- H
    fhat$gridtype <- "linear"
    fhat$gridded <- TRUE
    fhat$binned <- FALSE
    fhat$names <- parse.name(x)
    fhat$w <- rep(1, nrow(x))
    class(fhat) <- "kde"
                
    return(fhat)
}

dmvnorm.LB <- function(x, mu, Sigma, a0, a1)
{
    if (!is.matrix(x)) x <- as.matrix(x)
    if (missing(a0) | missing(a1))
    {
        ev <- -sweep(x, 2, mu, FUN="-") %*% diag(sqrt(1/diag(Sigma)))
        ev.list <- list(unique(ev[,1]), unique(ev[,2]))

        delta <- prod(unlist(lapply(lapply(ev.list, diff), getElement, 1)))
        d <- ncol(Sigma)
        evalK <- dmvnorm.mixt(x=ev, mus=rep(0,d), Sigmas=diag(d), props=1)
        m0 <- sum(evalK*delta)
        m1 <- apply(evalK*delta*ev, 2, sum)
        m2 <- apply(evalK*delta*rowKpow(ev, ev), 2, sum)
        M2 <- invvec(m2)
        M2inv <- chol2inv(chol(M2))
        a0 <- 1/drop(m0 - t(m1) %*% M2inv%*%m1)
        a1 <- -a0*M2inv%*%m1
    }    
    evalK.LB <- (a0 - drop(ev %*% a1))*dmvnorm.mixt(x=x, mus=mu, Sigmas=Sigma, props=1)

    return(evalK.LB)
}

######################################################################
## Boundary kernel estimator using beta bounday kernels (2nd form)
######################################################################

kde.beta.boundary <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, boundary.supp=1, verbose=FALSE)
{
    ## default values
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    if (missing(binned)) binned <- ksd$binned
    if (missing(gridsize)) gridsize <- ksd$gridsize
    if (missing(bgridsize)) bgridsize <- gridsize
    ## if (missing(gridsize)) gridsize <- default.gridsize(d)

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
    if (d==1 & missing(h)) 
        h <- hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n))
    if (d>1 & missing(H))
        H <- Hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n))
  
    ## compute exact (non-binned) estimator
    ## 1-dimensional    
    if (d==1)
    {
        if (missing(eval.points))
        {
            fhat <- kde.boundary.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned)
        }
        else
            stop("Not yet implemented.")  # fhat <- kde.points.1d(x=x, h=h, eval.points=eval.points, positive=positive, adj.positive=adj.positive, w=w)
    }
    ## multi-dimensional
    else
    {  
        if (is.data.frame(x)) x <- as.matrix(x)

        if (missing(eval.points))
        {
            if (d==2)
                fhat <- kde.boundary.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose)
            else if (d == 3)
                fhat <- kde.boundary.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, boundary.supp=boundary.supp, binned=binned, verbose=verbose) 
            else 
                stop("Need to specify eval.points for more than 3 dimensions")
        }
        else
            stop("Not yet implemented.")     
    }

    fhat$binned <- binned
    fhat$names <- parse.name(x)  ## add variable names
    fhat$w <- w
    class(fhat) <- "kde"
  
    ## compute prob contour levels
    if (compute.cont & missing(eval.points))
    fhat$cont <- contourLevels(fhat, cont=1:99, approx=approx.cont)

    return(fhat)
 }

kde.boundary.grid.1d <- function(x, h,gridsize, supp=3.7, xmin, xmax, gridtype, w, boundary.supp=0.5, binned=FALSE)
{
    if (missing(xmin)) xmin <- min(x)
    if (missing(xmax)) xmax <- max(x)
    if (missing(gridtype)) gridtype <- "linear"

    ## transform x into [0,1]
    x.star <- (x-xmin)/(xmax-xmin) 
    h.star <- h/(xmax-xmin)
    n <- length(x)

    gridtype1 <- match.arg(gridtype, c("linear", "sqrt")) 
    if (gridtype1=="linear")
        eval.x <- seq(0, 1, length=gridsize)
    else if (gridtype1=="sqrt")
    {
        eval.x.temp <- seq(0, 1, length=gridsize)
        eval.x <- sign(eval.x.temp) * eval.x.temp^2
    }
    gridtype.vec <- gridtype1

    ## indicator for closeness to boundary of [0,1]
    bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp) 

    n1 <- sum(!bound.ind)
    ## interior points - use normal kernel
    ## binned estimation only in the interior
    fhat.grid <- rep(0,length=gridsize)
    if (n1>0)
    {
        if (binned)
        {
            fhat.grid <- n1*kde(x=x.star[!bound.ind], h=h.star, xmin=0, xmax=1, binned=TRUE, bgridsize=gridsize)$estimate
        }
        else
        {
            fhat.grid <- n1*dnorm.mixt(x=eval.x, mus=x.star[!bound.ind], sigmas=rep(h.star, n1), props=w[!bound.ind]/n1)
        }
    }

    ## boundary points - use adjusted beta kernel
    hb.star <- 2*h.star
    for (i in 1:(n-n1))
        fhat.grid <- fhat.grid + dbeta.kernel2(x=x.star[bound.ind][i], eval.x=eval.x, h=hb.star)*w[bound.ind][i]
    fhat.grid <- fhat.grid/n

    ## backtransform
    eval.points <- (xmax-xmin)*eval.x + xmin 
    fhat.grid <- fhat.grid/(xmax-xmin)

    fhat <- list(x=x, eval.points=eval.points, estimate=fhat.grid, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE)
    class(fhat) <- "kde"

    return(fhat)
}

kde.boundary.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, boundary.supp=1, binned=FALSE, verbose=FALSE)
{
    n <- nrow(x)
    d <- ncol(x)

    if (missing(xmin)) xmin <- apply(x, 2, min) ##- h*supp
    if (missing(xmax)) xmax <- apply(x, 2, max) ## + h*supp
    if (missing(gridtype)) gridtype <- rep("linear", d)

    ## transform x into [0,1]^d
    x.star <- x
    for (j in 1:d) x.star[,j] <- (x[,j]-xmin[j])/(xmax[j]-xmin[j])
    H.star <- diag(1/(xmax-xmin)) %*% H %*% diag(1/(xmax-xmin))
    h.star <- sqrt(diag(H.star))

    ## initialise grid 
    if (is.null(gridx))
        gridx <- make.grid.ks(x.star, matrix.sqrt(H.star), tol=supp, gridsize=gridsize, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype) 
    suppx <- make.supp(x.star, matrix.sqrt(H.star), tol=supp)

    if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)    
    fhat.grid <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))

    ## indicator for closeness to boundary of [0,1]^d
    bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp)
    n1 <- sum(!bound.ind)

    if (verbose) pb <- txtProgressBar()
    if (binned)
    {
        ## interior points - use normal kernel
        fhat.grid <- n1*kde(x=x.star[!bound.ind,], H=H.star, xmin=rep(0,d), xmax=rep(1,d), binned=TRUE, bgridsize=gridsize)$estimate

        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)
            if (bound.ind[i])
            {
                ## compute evaluation points 
                eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
                eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
                eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
                eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
                eval.x.len <- length(eval.x)

                ## convert bandwidth from normal kernel to beta kernel scale
                ## for boundary points 
                hb.star <- 2*h.star
                fhat <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star)

                ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
                for (j in 1:length(eval.y))
                    fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
                        fhat.grid[eval.x.ind, eval.y.ind[j]] + 
                        w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
            }
        }
    }
    else
    {  
        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)

            ## compute evaluation points 
            eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
            eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
            eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
            eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
            eval.x.len <- length(eval.x)

            ## interior points - use normal kernel
            if (!bound.ind[i])
            {
            eval.pts <- expand.grid(eval.x, eval.y)
            fhat <- dmvnorm(eval.pts, x.star[i,], H.star)
            }
            else
            {
            ## convert bandwidth from normal kernel to beta kernel scale
            ## for boundary points 
            hb.star <- 2*h.star
            fhat <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star)
        }

        ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        for (j in 1:length(eval.y))
            fhat.grid[eval.x.ind, eval.y.ind[j]] <- 
                fhat.grid[eval.x.ind, eval.y.ind[j]] + 
                w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
        }
    }
    fhat.grid <- fhat.grid/n
    if (verbose) close(pb)
  
    ## back-transform
    gridx1 <- list((xmax[1]-xmin[1])*gridx[[1]] + xmin[1], (xmax[2]-xmin[2])*gridx[[2]] + xmin[2]) 
    fhat.grid <- fhat.grid/prod(xmax-xmin)

    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, boundary=bound.ind)

    return(fhat.list)
}

kde.boundary.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, boundary.supp=0.5, verbose=FALSE, binned=FALSE)
{
    n <- nrow(x)
    d <- ncol(x)

    if (missing(xmin)) xmin <- apply(x, 2, min) 
    if (missing(xmax)) xmax <- apply(x, 2, max) 
    if (missing(gridtype)) gridtype <- rep("linear", d)

    ## transform x into [0,1]^d
    x.star <- x
    for (j in 1:d) x.star[,j] <- (x[,j]-xmin[j])/(xmax[j]-xmin[j])
    H.star <- diag(1/(xmax-xmin)) %*% H %*% diag(1/(xmax-xmin))
    h.star <- sqrt(diag(H.star)) 

    ## initialise grid 
    if (is.null(gridx))
        gridx <- make.grid.ks(x.star, matrix.sqrt(H.star), tol=supp, gridsize=gridsize, xmin=rep(0,d), xmax=rep(1,d), gridtype=gridtype) 
    suppx <- make.supp(x.star, matrix.sqrt(H.star), tol=supp)
    if (is.null(grid.pts)) grid.pts <- find.gridpts(gridx, suppx)
    fhat.grid <- array(0, dim=c(length(gridx[[1]]), length(gridx[[2]]), length(gridx[[3]])))

    ## indicator for closeness to boundary of [0,1]^d
    bound.ind <- boundary.ind(x=x.star, h=h.star, boundary.supp=boundary.supp)
    n1 <- sum(!bound.ind)

    if (verbose) pb <- txtProgressBar() 
    if (binned)
    {
        ## interior points - use normal kernel
        fhat.grid <- n1*kde(x=x.star[!bound.ind,], H=H.star, xmin=rep(0,d), xmax=rep(1,d), binned=TRUE, bgridsize=gridsize)$estimate

        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)

            if (bound.ind[i])
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

                ## convert bandwidth from normal kernel to beta kernel scale
                ## for boundary points 
                hb.star <- 2*h.star
                fhat.xy <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star[1:2])

                ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
                for (k in 1:length(eval.z))
                {
                    fhat <- w[i]*cbind(fhat.xy, dbeta.kernel2(x=x.star[i,3], eval.x=eval.z[k], h=hb.star[3]))
                    for (j in 1:length(eval.y))
                    fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                      fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                    fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
                }
            }
        }
    }
    else
    {  
        for (i in 1:n)
        {
            if (verbose) setTxtProgressBar(pb, i/n)

            ## compute evaluation points
            eval.x <- gridx[[1]][grid.pts$xmin[i,1]:grid.pts$xmax[i,1]]
            eval.y <- gridx[[2]][grid.pts$xmin[i,2]:grid.pts$xmax[i,2]]
            eval.z <- gridx[[3]][grid.pts$xmin[i,3]:grid.pts$xmax[i,3]]
            eval.x.ind <- c(grid.pts$xmin[i,1]:grid.pts$xmax[i,1])
            eval.y.ind <- c(grid.pts$xmin[i,2]:grid.pts$xmax[i,2])
            eval.z.ind <- c(grid.pts$xmin[i,3]:grid.pts$xmax[i,3])
            eval.x.len <- length(eval.x)
            eval.pts <- expand.grid(eval.x, eval.y)

            ## interior points - use normal kernel
            if (!bound.ind[i])
            {
                ## place vector of density estimate values `fhat' onto grid 'fhat.grid' 
        
                for (k in 1:length(eval.z))
                {
                    fhat <- w[i]*dmvnorm(cbind(eval.pts, eval.z[k]), x[i,], H)
                    for (j in 1:length(eval.y))
                        fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                            fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                            fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
                }
            }
            else
            {
                ## convert bandwidth from normal kernel to beta kernel scale
                hb.star <- 2*h.star
                fhat.xy <- dmvbeta.prod.kernel2(x=x.star[i,], eval.x=list(eval.x, eval.y), hs=hb.star[1:2])
                    
                for (k in 1:length(eval.z))
                {
                    fhat <- w[i]*cbind(fhat.xy, dbeta.kernel2(x=x.star[i,3], eval.x=eval.z[k], h=hb.star[3]))
                    for (j in 1:length(eval.y))
                        fhat.grid[eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                            fhat.grid[eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                            fhat[((j-1) * eval.x.len + 1):(j * eval.x.len)]
                }
            }
        }
    }
    fhat.grid <- fhat.grid/n
    if (verbose) close(pb)
  
    ## back-transform
    gridx1 <- list((xmax[1]-xmin[1])*gridx[[1]] + xmin[1], (xmax[2]-xmin[2])*gridx[[2]] + xmin[2], (xmax[3]-xmin[3])*gridx[[3]] + xmin[3])
    fhat.grid <- fhat.grid/prod(xmax-xmin)

    fhat.list <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, boundary=bound.ind)
    
    return(fhat.list)
}

## indicator function for boundary region of [0,1] i.e. [0,h] + [1-h, h]

boundary.ind <- function(x, h, xmin, xmax, boundary.supp=1)
{
    if (is.vector(x)) { x <- matrix(x, ncol=1) }
    n <- nrow(x)
    d <- ncol(x)

    ## indicator for closeness to boundary of [0,1]^d
    bound.ind <- matrix(NA, nrow=n, ncol=d)
    for (j in 1:d) bound.ind[,j] <- (abs(x[,j]) <= boundary.supp*h[j]) | (abs(1-x[,j]) <= boundary.supp*h[j])
    bound.ind <- apply(bound.ind, 1, any)

    return(bound.ind)
}

######################################################################
## Bivariate beta boundary KDE
######################################################################

## modified boundary beta kernel - first form (Chen, 1999)
dbeta.kernel <- function(x, eval.x, h)
{
    return (dbeta(x=eval.x, shape1=x/h^2+1, shape2=(1-x)/h^2+1))
}

## modified boundary beta kernel - second form (Chen, 1999)
dbeta.kernel2 <- function(x, eval.x, h)
{
    rhox <- function(y, hy) { if (y==0) return (1) else return(2*hy^4 + 5/2 - sqrt(4*hy^8 + 6*hy^4 + 9/4 - y^2 -y/hy^2)) }
    ind <- cut(x, c(0, 2*h^2, 1-2*h^2, 1), labels=FALSE, include.lowest=TRUE)
    dbk <- rep(0, length(eval.x))

    if (ind==1) { shape1 <- rhox(x,hy=h); shape2 <- (1-x)/h^2 }
    else if (ind==2) { shape1 <- x/h^2;   shape2 <- (1-x)/h^2 }
    else if (ind==3) { shape1 <- x/h^2;   shape2 <- rhox(1-x, hy=h) }

return(dbeta(eval.x, shape1=shape1, shape2=shape2))
}

## modified multivariate boundary beta product kernel

dmvbeta.prod.kernel2 <- function(x, eval.x, hs)
{
    d <- length(hs)
    db <- vector("list", d)
    for (i in 1:d)
    {   
        db[[i]] <- 0
        db[[i]] <- dbeta.kernel2(x=x[i], eval.x=eval.x[[i]], h=hs[i])
    }

    db <- expand.grid(db)
    db <- apply(db, 1, prod)
    
    return(db)
}

## modified multivariate boundary beta spherically symmetric kernel

dmvbeta.symm.kernel2 <- function(x, eval.x, H)
{
    d <- ncol(H)
    eval.y <- sqrt(apply(eval.x^2, 1, sum))/sqrt(d)
    y <- sqrt(sum(x^2))/sqrt(d)
  
    return(dbeta.kernel2(x=y, eval.x=eval.y, h=sqrt(tr(H)))/d)
}

rbeta.kernel2 <- function(x, n, h)
{
    rhox <- function(y, hy) { if (y==0) return (1) else return(2*hy^4 + 5/2 - sqrt(4*hy^8 + 6*hy^4 + 9/4 - y^2 -y/hy^2)) }
    
    ind <- cut(x, c(0, 2*h^2, 1-2*h^2, 1), labels=FALSE, include.lowest=TRUE)

    if (ind==1) { shape1 <- rhox(x,hy=h); shape2 <- (1-x)/h^2 }
    else if (ind==2) { shape1 <- x/h^2;   shape2 <- (1-x)/h^2 }
    else if (ind==3) { shape1 <- x/h^2;   shape2 <- rhox(1-x, hy=h) }
    
    return(rbeta(n=n, shape1=shape1, shape2=shape2))
}

dmvbeta.prod.kernel2.2d <- function(x, hs, xmin=c(0,0), xmax=c(1,1), ...)
{
    x.star <- (x - xmin)/(xmax - xmin)
    hs.star <- hs/(xmax - xmin)
    x1 <- seq(0,1, length=151)
    x1[1] <- 1e-9; x1[151] <- 1-1e-9
    eval.points <- list(x1, x1)
    fhat <- list()
    x <- cbind(rbeta.kernel2(x=x.star[1], n=1000, h=hs.star[1]), rbeta.kernel2(x=x.star[2], n=1000, h=hs.star[2]))
   
    fhat$eval.points <- list(seq(0,1, length=151), seq(0,1, length=151))
    
    fhat$estimate <- matrix(dmvbeta.prod.kernel2(x=x.star, eval.x=eval.points, hs=hs.star, ...), nrow=length(eval.points[[1]]))

    fhat$x <- sweep(sweep(x, 2, xmin, FUN="+"), 2, xmax-xmin, FUN="*") 
    fhat$eval.points[[1]] <- xmin[1] + fhat$eval.points[[1]]*(xmax[1]-xmin[1])
    fhat$eval.points[[2]] <- xmin[2] + fhat$eval.points[[2]]*(xmax[2]-xmin[2])
    fhat$H <- diag(2)
    fhat$gridtype <- "linear"
    fhat$gridded <- TRUE
    fhat$binned <- FALSE
    fhat$names <- parse.name(x)
    fhat$w <- rep(1, nrow(x))
    class(fhat) <- "kde"
                
    return(fhat)
}

##########################################################################
## Truncate unbounded KDE to polygon boundary
##########################################################################

kde.truncate <- function(fhat, boundary)
{
    ## reallocate any probability mass outside of map boundary regions
    ## to interior regions
    truncate.ind <- array(mgcv::in.out(boundary, as.matrix(expand.grid(fhat$eval.points[[1]], fhat$eval.points[[2]]))), dim=dim(fhat$estimate))
    fhat.trunc <- fhat
    fhat.trunc$estimate <- fhat.trunc$estimate*truncate.ind
    fhat.sum <- contourProbs(fhat, abs.cont=0)
    fhat.trunc.sum <- contourProbs(fhat.trunc, abs.cont=0)
    fhat.trunc$estimate <- fhat.trunc$estimate*fhat.sum/fhat.trunc.sum
    
    fhat.trunc.temp <- fhat.trunc
    fhat.trunc.temp$x <- fhat.trunc.temp$x[mgcv::in.out(boundary, fhat.trunc$x),]
    fhat.trunc$cont <- contourLevels(fhat.trunc.temp, cont=1:99, approx=TRUE)
    
    return(fhat.trunc)
}

##########################################################################
## Truncate unbounded KDDE to polygon boundary
##########################################################################

kdde.truncate <- function(fhat, boundary)
{
    ## reallocate any probability mass outside of map boundary regions
    ## to interior regions
    fhat.trunc <- fhat
    fhat0 <- kde(x=fhat$x, H=fhat$H)
    fhat0.sum <- contourProbs(fhat0, abs.cont=0)
    fhat0.trunc <- kde.truncate(fhat=fhat0, boundary=boundary)
    fhat0.trunc.sum <- contourProbs(fhat0.trunc, abs.cont=0)
    
    for (i in 1:length(fhat$estimate))
    {    
        truncate.ind <- array(mgcv::in.out(boundary, as.matrix(expand.grid(fhat$eval.points[[1]], fhat$eval.points[[2]]))), dim=dim(fhat$estimate[[i]]))
        fhat.trunc$estimate[[i]] <- fhat.trunc$estimate[[i]]*truncate.ind
        fhat.trunc$estimate[[i]] <- fhat.trunc$estimate[[i]]*fhat0.sum/fhat0.trunc.sum
    }
    fhat.trunc.temp <- fhat.trunc
    fhat.trunc.temp$x <- fhat.trunc.temp$x[mgcv::in.out(boundary, fhat.trunc$x),]
    fhat.trunc$cont <- contourLevels(fhat.trunc.temp, cont=1:99, approx=TRUE)
    
    return(fhat.trunc)
}
