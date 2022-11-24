###############################################################################
## Multivariate kernel density derivative estimate 
###############################################################################

kdde <- function(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, positive=FALSE, adj.positive, w, deriv.vec=TRUE, verbose=FALSE)
{
    ## default values 
    r <- deriv.order
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize

    ## clip data to xmin,xmax grid for binned estimation
    grid.clip <- binned 
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        if (positive & missing(xmin)) { xmin <- rep(0,d) }
        xt <- truncate.grid(x=x, y=w, xmin=xmin, xmax=xmax)
        x <- xt$x; w <- xt$y; n <- length(w)
    }
  
    ## default bandwidths
    if (d==1 & missing(h)) 
    {
        if (positive) x1 <- log(x) else x1 <- x
        h <- hpi(x=x1, nstage=2, binned=default.bflag(d=d, n=n), deriv.order=r)
    }
    if (d>1 & missing(H))
    {
        if (positive) x1 <- log(x) else x1 <- x
        if ((r>0) & (d>2)) nstage <- 1 else nstage <- 2
        H <- Hpi(x=x1, nstage=nstage, binned=default.bflag(d=d, n=n), deriv.order=r, verbose=verbose)
    }
  
    ## compute binned estimator
    if (binned)
    {
        if (positive & is.vector(x))
        {
            y <- log(x)
            fhat <- kdde.binned(x=y, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w)
            fhat$estimate <- fhat$estimate/exp(fhat$eval.points)
            fhat$eval.points <- exp(fhat$eval.points)
            fhat$x <- x
        }
        else
            fhat <- kdde.binned(x=x, H=H, h=h, deriv.order=r, bgridsize=bgridsize, xmin=xmin, xmax=xmax, w=w, deriv.vec=deriv.vec, verbose=verbose)

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
            if (!missing(H) & !missing(h)) stop("Both H and h are both specified")
          
            if (missing(h)) h <- sqrt(H)

            if (missing(eval.points))
                fhat <- kdde.grid.1d(x=x, h=h, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r)
            else
                fhat <- kdde.points.1d(x=x, h=h, eval.points=eval.points, w=w, deriv.order=r)
        }
        ## multi-dimensional
        else
        {   
            if (is.data.frame(x)) x <- as.matrix(x)
          
            if (missing(eval.points))
            {
                if (d==2) 
                    fhat <- kdde.grid.2d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r, deriv.vec=deriv.vec, verbose=verbose)
                else if (d==3)
                    fhat <- kdde.grid.3d(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, deriv.order=r, deriv.vec=deriv.vec, verbose=verbose) 
                else 
                    stop("Need to specify eval.points for more than 3 dimensions")
            }
            else
                fhat <- kdde.points(x=x, H=H, eval.points=eval.points, w=w, deriv.order=r, deriv.vec=deriv.vec)
        }
    }

    fhat$binned <- binned
    fhat$names <- parse.name(x)
    fhat$type <- "kdde"
    class(fhat) <- "kdde"
    
    return(fhat)
}

###############################################################################
## Multivariate binned kernel density derivative estimate
###############################################################################

kdde.binned <- function(x, H, h, deriv.order, bgridsize, xmin, xmax, bin.par, w, deriv.vec=TRUE, deriv.index, verbose=FALSE)
{
    r <- deriv.order
    if (length(r)>1) stop("deriv.order should be a non-negative integer")

    ## linear binning
    if (missing(bin.par))
    {
        if (is.vector(x)) { d  <- 1 }
        else { d <- ncol(x) }

        if (d==1)
            if (missing(H)) { H <- as.matrix(h^2) } 
            else { h <- sqrt(H); H <- as.matrix(H) }

        if (d==1) Hd <- H else Hd <- diag(diag(H)) 
        bin.par <- binning(x=x, H=Hd, h=h, bgridsize=bgridsize, xmin=xmin, xmax=xmax, supp=3.7+max(r), w=w)
    }
    else
    {
        if (!is.list(bin.par$eval.points)) { d <- 1; bgridsize <- length(bin.par$eval.points) }
        else  { d <- length(bin.par$eval.points); bgridsize <- sapply(bin.par$eval.points, length) } 

        w <- bin.par$w
        if (d==1)
            if (missing(H)) H <- as.matrix(h^2)
            else { h <- sqrt(H); H <- as.matrix(H) }
    }
     
    if (d==1)
    {
        fhat <- kdde.binned.1d(h=h, deriv.order=r, bin.par=bin.par)
        eval.points <- fhat$eval.points
        est <- fhat$estimate
    }
    else
    {
        ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE, deriv.vec=deriv.vec)
        fhat.grid <- kdde.binned.nd(H=H, deriv.order=r, bin.par=bin.par, verbose=verbose, deriv.vec=deriv.vec)
    }

    if (missing(x)) x <- NULL
  
    if (d==1)
    {
        if (r==0) fhat <- list(x=x, eval.points=unlist(eval.points), estimate=est, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w)
        else
        fhat <- list(x=x, eval.points=unlist(eval.points), estimate=est, h=h, H=h^2, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=r)
      }
  else
  {
        if (r==0)
            fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate[[1]], H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w)
        else
            fhat <- list(x=x, eval.points=fhat.grid$eval.points, estimate=fhat.grid$estimate, H=H, gridtype="linear", gridded=TRUE, binned=TRUE, names=NULL, w=w, deriv.order=r, deriv.ind=ind.mat)
  }
  class(fhat) <- "kdde"
  
  return(fhat)
}

kdde.binned.1d <- function(h, deriv.order, bin.par)
{
    r <- deriv.order
    n <- sum(bin.par$counts)
    a <- min(bin.par$eval.points)
    b <- max(bin.par$eval.points)
    M <- length(bin.par$eval.points)
    L <- min(ceiling((4+r)*h*(M-1)/(b-a)), M-1)
    delta <- (b-a)/(M-1)
    N <- 2*L-1
    grid1 <- seq(-(L-1), L-1)

    keval <- dnorm.deriv(x=delta*grid1, mu=0, sigma=h, deriv.order=r)/n
    est <- symconv.1d(keval, bin.par$counts)

    return(list(eval.points=bin.par$eval.points, estimate=est))
}

kdde.binned.nd <- function(H, deriv.order, bin.par, verbose=FALSE, deriv.vec=TRUE)
{
    d <- ncol(H)
    r <- deriv.order
    n <- sum(bin.par$counts)
    a <- sapply(bin.par$eval.points, min)
    b <- sapply(bin.par$eval.points, max)
    M <- sapply(bin.par$eval.points, length)
    L <- pmin(ceiling((4+r)*max(sqrt(abs(diag(H))))*(M-1)/(b-a)), M-1)
    delta <- (b-a)/(M-1)

    if (min(L)<=0) warning(paste("Binning grid too coarse for current (small) bandwidth: consider increasing grid size for dimensions", toString(which(pmin(L)<=1))))
    N <- 2*L-1
    if(d==2)
    {
        grid1 <- seq(-(L[1]-1), L[1]-1)
        grid2 <- seq(-(L[2]-1), L[2]-1)
        xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2)
    }
    else if (d==3)
    {
        grid1 <- seq(-(L[1]-1), L[1]-1)
        grid2 <- seq(-(L[2]-1), L[2]-1)
        grid3 <- seq(-(L[3]-1), L[3]-1)
        xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2, delta[3]*grid3)
    }
    else if (d==4)
    {
        grid1 <- seq(-(L[1]-1), L[1]-1)
        grid2 <- seq(-(L[2]-1), L[2]-1)
        grid3 <- seq(-(L[3]-1), L[3]-1)
        grid4 <- seq(-(L[4]-1), L[4]-1)
        xgrid <- expand.grid(delta[1]*grid1, delta[2]*grid2, delta[3]*grid3, delta[4]*grid4)
    }
    deriv.index <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, only.index=TRUE, deriv.vec=TRUE) 
    deriv.index.minimal <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, only.index=TRUE, deriv.vec=FALSE)

    if (verbose) pb <- txtProgressBar()

    if (r==0)
    {
        n.seq <- block.indices(1, nrow(xgrid), d=d, r=r, diff=FALSE)
        est.list <- vector(1, mode="list")
        est.list[[1]] <- array(0, dim=dim(bin.par$counts))
    }
    else if (r>0)
    {      
          n.deriv <- nrow(deriv.index)
          n.deriv.minimal <- nrow(deriv.index.minimal)
          if (deriv.vec) n.est.list <- n.deriv else n.est.list <- n.deriv.minimal
          est.list <- vector(n.est.list, mode="list")
          
          for (j in 1:n.est.list) est.list[[j]] <- array(0, dim=dim(bin.par$counts))
          if (d^r >= 3^7) n.seq <- block.indices(1, nrow(xgrid), d=d, r=r, diff=FALSE, block.limit=1e4)
          else n.seq <- block.indices(1, nrow(xgrid), d=d, r=r, diff=FALSE, block.limit=1e5)
    }

    for (i in 1:(length(n.seq)-1))
    {
        if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
        keval <- dmvnorm.deriv(x=xgrid[n.seq[i]:(n.seq[i+1]-1),], mu=rep(0,d), Sigma=H, deriv.order=r, add.index=TRUE, deriv.vec=FALSE)$deriv/n

        if (r==0) keval <- as.matrix(keval, ncol=1)
        else if (is.vector(keval)) keval <- as.matrix(t(keval), nrow=1)
        est <- list()

        ## loop over only unique partial derivative indices
        nderiv <- nrow(deriv.index.minimal)
        if (!(is.null(nderiv)))
            for (s in 1:nderiv)
            {
                if (deriv.vec) deriv.rep.index <- which.mat(deriv.index.minimal[s,], deriv.index)
                else deriv.rep.index <- s
                kevals <- array(keval[,s], dim=N)
                est.temp <- symconv.nd(kevals, bin.par$counts, d=d)
                for (s2 in 1:length(deriv.rep.index)) est[[deriv.rep.index[s2]]] <- est.temp
            }
        else
        {
            kevals <- lapply(as.data.frame(keval), function(x) { array(x, dim=N) })
            est <- lapply(kevals, function(x) { symconv.nd(x, bin.par$counts, d=d) })
        }
        if (r==0) est.list[[1]] <- est.list[[1]] + est[[1]]
        else if (r>0) for (j in 1:n.est.list) est.list[[j]] <- est.list[[j]] + est[[j]]  
    }
    if (verbose) close(pb)
    fhatr <- list(eval.points=bin.par$eval.points, estimate=est.list, deriv.order=r)

    return(fhatr)
}

#############################################################################
## Univariate kernel density derivative estimate on a grid
#############################################################################

kdde.grid.1d <- function(x, h, gridsize, supp=3.7, positive=FALSE, adj.positive, xmin, xmax, gridtype, w, deriv.order=0)
{
    r <- deriv.order
    if (r==0)
        fhatr <- kde(x=x, h=h, gridsize=gridsize, supp=supp, positive=positive, adj.positive=adj.positive, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w)
    else
    {  
        if (missing(xmin)) xmin <- min(x) - h*supp
        if (missing(xmax)) xmax <- max(x) + h*supp
        if (missing(gridtype)) gridtype <- "linear"

        y <- x
        gridtype1 <- match.arg(gridtype, c("linear", "sqrt"))
        if (gridtype1=="linear")
            gridy <- seq(xmin, xmax, length=gridsize)
        else if (gridtype1=="sqrt")
        {
            gridy.temp <- seq(sign(xmin)*sqrt(abs(xmin)), sign(xmax)*sqrt(abs(xmax)), length=gridsize)
            gridy <- sign(gridy.temp) * gridy.temp^2
        }
        gridtype.vec <- gridtype1

        n <- length(y)
        est <- dnorm.deriv.mixt(x=gridy, mus=y, sigmas=rep(h, n), props=w/n, deriv.order=r)
        fhatr <- list(x=y, eval.points=gridy, estimate=est, h=h, H=h^2, gridtype=gridtype.vec, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=deriv.order)
        class(fhatr) <- "kdde"
    }

    return(fhatr)
}

##############################################################################
## Bivariate kernel density derivative estimate on a grid
## Computes all mixed partial derivatives for a given deriv.order
##############################################################################

kdde.grid.2d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0, deriv.vec=TRUE, verbose=FALSE)
{
    d <- 2
    r <- deriv.order
    if (r==0)
    fhatr <- kde(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
    else
    {  
        ## initialise grid 
        n <- nrow(x)
        if (is.null(gridx))
          gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 
        
        suppx <- make.supp(x, matrix.sqrt(H), tol=supp)
        
        if (is.null(grid.pts))
        grid.pts <- find.gridpts(gridx, suppx)    

        nderiv <- d^r
        fhat.grid <- list()
        for (k in 1:nderiv)
            fhat.grid[[k]] <- matrix(0, nrow=length(gridx[[1]]), ncol=length(gridx[[2]]))
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

            ## Create list of matrices for different partial derivatives
            fhat <- dmvnorm.deriv(x=eval.pts, mu=x[i,], Sigma=H, deriv.order=r)

            ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
            for (k in 1:nderiv)
            for (j in 1:length(eval.y))
                fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] <- fhat.grid[[k]][eval.x.ind, eval.y.ind[j]] + w[i]*fhat[((j-1) * eval.x.len + 1):(j * eval.x.len),k]
            if (verbose) setTxtProgressBar(pb, i/n)
        }
        if (verbose) close(pb)
        
        for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
        gridx1 <- list(gridx[[1]], gridx[[2]]) 

        ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE)

        if (!deriv.vec)
        {
            fhat.grid.vech <- list()
            deriv.ind <- unique(ind.mat)

            for (i in 1:nrow(deriv.ind))
            {
                which.deriv <- which.mat(deriv.ind[i,], ind.mat)[1]
                fhat.grid.vech[[i]] <- fhat.grid[[which.deriv]]
            }
            ind.mat <- deriv.ind
            fhat.grid <- fhat.grid.vech
        }

        fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
      }

      return(fhatr)
}

##############################################################################
## Trivariate kernel density derivative estimate on a grid
## Computes all mixed partial derivatives for a given deriv.order
##############################################################################

kdde.grid.3d <- function(x, H, gridsize, supp, gridx=NULL, grid.pts=NULL, xmin, xmax, gridtype, w, deriv.order=0, deriv.vec=TRUE, verbose=FALSE)
{
    d <- 3
    r <- deriv.order
    if (r==0)
        fhatr <- kde(x=x, H=H, gridsize=gridsize, supp=supp, xmin=xmin, xmax=xmax, gridtype=gridtype, w=w, verbose=verbose)
    else
    {
        ## initialise grid 
        n <- nrow(x)
        if (is.null(gridx))
            gridx <- make.grid.ks(x, matrix.sqrt(H), tol=supp, gridsize=gridsize, xmin=xmin, xmax=xmax, gridtype=gridtype) 

        suppx <- make.supp(x, matrix.sqrt(H), tol=supp)

        if (is.null(grid.pts))
        grid.pts <- find.gridpts(gridx, suppx)    

        nderiv <- d^r
        fhat.grid <- list()
        for (k in 1:nderiv)
            fhat.grid[[k]] <- array(0, dim=gridsize)
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

            ## create list of matrices for different partial derivatives
            ## place vector of density estimate values `fhat' onto grid 'fhat.grid'
            for (ell in 1:nderiv)
                for (k in 1:length(eval.z))
                {
                    fhat <- w[i]*dmvnorm.deriv(cbind(eval.pts, eval.z[k]), x[i,], H, deriv.order=r)
                    for (j in 1:length(eval.y))
                        fhat.grid[[ell]][eval.x.ind,eval.y.ind[j], eval.z.ind[k]] <- 
                            fhat.grid[[ell]][eval.x.ind, eval.y.ind[j], eval.z.ind[k]] + 
                            fhat[((j-1) * eval.x.len + 1):(j * eval.x.len), ell]
                }
          
            if (verbose) setTxtProgressBar(pb, i/n)
        }
        if (verbose) close(pb)

        for (k in 1:nderiv) fhat.grid[[k]] <- fhat.grid[[k]]/n
        gridx1 <- list(gridx[[1]], gridx[[2]], gridx[[3]])

        ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=H, deriv.order=r, only.index=TRUE)

        if (!deriv.vec)
        {
            fhat.grid.vech <- list()
            deriv.ind <- unique(ind.mat)

            for (i in 1:nrow(deriv.ind))
            {
                which.deriv <- which.mat(deriv.ind[i,], ind.mat)[1]
                fhat.grid.vech[[i]] <- fhat.grid[[which.deriv]]
            }
            ind.mat <- deriv.ind
            fhat.grid <- fhat.grid.vech
        }

        fhatr <- list(x=x, eval.points=gridx1, estimate=fhat.grid, H=H, gridtype=gridx$gridtype, gridded=TRUE, binned=FALSE, names=NULL, w=w, deriv.order=deriv.order, deriv.ind=ind.mat)
    }

    return(fhatr)
}

#############################################################################
## Multivariate kernel density estimate using normal kernels,
## evaluated at each sample point
#############################################################################

kdde.points.1d <- function(x, h, eval.points, w, deriv.order=0) 
{
    r <- deriv.order
    n <- length(x)
    fhat <- dnorm.deriv.mixt(x=eval.points, mus=x, sigmas=rep(h,n), props=w/n, deriv.order=r)

    return(list(x=x, eval.points=eval.points, estimate=fhat, h=h, H=h^2, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=r))
}

kdde.points <- function(x, H, eval.points, w, deriv.order=0, deriv.vec=TRUE) 
{
    n <- nrow(x)
    Hs <- replicate(n, H, simplify=FALSE) 
    Hs <- do.call(rbind, Hs)
    r <- deriv.order
    fhat <- dmvnorm.deriv.mixt(x=eval.points, mus=x, Sigmas=Hs, props=w/n, deriv.order=r, deriv.vec=deriv.vec, add.index=TRUE)

    return(list(x=x, eval.points=eval.points, estimate=fhat$deriv, H=H, gridded=FALSE, binned=FALSE, names=NULL, w=w, deriv.order=r, deriv.ind=fhat$deriv.ind))
}

#############################################################################
## S3 methods for KDDE objects
#############################################################################

## plot method

plot.kdde <- function(x, ...)
{
    fhat <- x

    if (is.null(fhat$deriv.order))
    {
        class(fhat) <- "kde"
        plot(fhat, ...)
    }
    else
    {  
        if (is.vector(fhat$x))
        {
            plotkdde.1d(fhat, ...)
            invisible()
        }
        else
        {
            d <- ncol(fhat$x)
            if (d==2) 
            {
                opr <- options()$preferRaster; if (!is.null(opr)) if (!opr) options("preferRaster"=TRUE)
                plotret <- plotkdde.2d(fhat, ...)
                if (!is.null(opr)) options("preferRaster"=opr)
                invisible(plotret)
            }
            else if (d==3)
            {
                plotkdde.3d(fhat, ...)
                invisible()
            }
            else 
                stop ("Plot function only available for 1, 2 or 3-d data")
        }
    }
}

plotkdde.1d <- function(fhat, xlab, ylab="Density derivative function", cont=50, abs.cont, ...)
{
    if (missing(xlab)) xlab <- fhat$names[1]
    
    if (missing(abs.cont))
    {
        abs.cont <- as.matrix(contourLevels(fhat, approx=TRUE, cont=cont), ncol = length(cont))
        abs.cont <- c(abs.cont[1, ], rev(abs.cont[2, ]))
    }
    class(fhat) <- "kde"

    plot(fhat, xlab=xlab, ylab=ylab, abs.cont=abs.cont, ...)
}

plotkdde.2d <- function(fhat, which.deriv.ind=1, cont=c(25,50,75), abs.cont, display="slice", xlab, ylab, zlab="Density derivative function", col, col.fun, alpha=1, kdde.flag=TRUE, thin=3, transf=1, neg.grad=FALSE, ...)
{
    disp <- match.arg(display, c("slice", "persp", "image", "filled.contour", "filled.contour2", "quiver"))
    if (disp=="filled.contour2") disp <- "filled.contour"
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Blue-Red", alpha=alpha) }
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (disp=="slice" | disp=="filled.contour")
    {
        if (missing(abs.cont))
        {
            abs.cont <- as.matrix(contourLevels(fhat, approx=TRUE, cont=cont, which.deriv.ind=which.deriv.ind), ncol=length(cont))
            abs.cont <- c(abs.cont[1,], rev(abs.cont[2,]))
        }
        if (missing(col))
        {
            if (disp=="slice") col <- col.fun(n=length(abs.cont)) 
            else if (disp=="filled.contour") 
            {
                col <- col.fun(n=length(abs.cont)+1)
                col[median(1:length(col))] <- "transparent"
            }
        }
    } 

    if (disp=="quiver")
    {
        if (fhat$deriv.order==1)
            plotquiver(fhat=fhat, thin=thin, transf=transf, neg.grad=neg.grad, col=col, xlab=xlab, ylab=ylab, alpha=alpha, ...)
        else warning("Quiver plot requires gradient estimate.")
    }
    else
    {
        fhat.temp <- fhat
        fhat.temp$deriv.ind <- fhat.temp$deriv.ind[which.deriv.ind,]
        fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
        fhat <- fhat.temp
        class(fhat) <- "kde"
      
        if (disp=="persp") plot(fhat, display=display, abs.cont=abs.cont, xlab=xlab, ylab=ylab, zlab=zlab, col.fun=col.fun, kdde.flag=kdde.flag, col=col, thin=thin, alpha=alpha, ...) 
        else plot(fhat, display=display, abs.cont=abs.cont, xlab=xlab, ylab=ylab, zlab=zlab, col.fun=col.fun, kdde.flag=kdde.flag, col=col, alpha=alpha, ...) 
    }
}

plotkdde.3d <- function(fhat, which.deriv.ind=1, display="plot3D", cont=c(25,50,75), abs.cont, colors, col, col.fun, xlab, ylab, zlab, ...)
{
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Blue-Red") }
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (missing(zlab)) ylab <- fhat$names[3]

    if (missing(abs.cont))
    {
        abs.cont <- as.matrix(contourLevels(fhat, approx=TRUE, cont=cont, which.deriv.ind=which.deriv.ind), ncol=length(cont))
        abs.cont <- c(abs.cont[1,], rev(abs.cont[2,]))
    }

    fhat.temp <- fhat 
    fhat.temp$deriv.ind <- fhat.temp$deriv.ind[which.deriv.ind,]
    fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
    fhat <- fhat.temp
    class(fhat) <- "kde"

    if (missing(col))
    {
        col <- col.fun(n=length(abs.cont)+1)
        nc <- length(col)
        col <- col[-median(1:nc)] 
    }

    colors <- col
    plot(fhat, display=display, abs.cont=abs.cont, colors=col, col=col, xlab=xlab, ylab=ylab, zlab=zlab, ...) 
}

######################################################################
## Quiver plot
######################################################################

plotquiver <- function(fhat, thin=5, transf=1, neg.grad=FALSE, xlab, ylab, col, add=FALSE, scale, length=0.1, alpha=1, ...)
{
    if (!requireNamespace("pracma", quietly=TRUE)) stop("Install the pracma package as it is required.", call.=FALSE)
    if (missing(col)) col <- 1
    col <- transparency.col(col, alpha=alpha)
    
    ev <- fhat$eval.points
    est <- fhat$estimate

    if (transf!=0)
    {
        est[[1]] <- sign(est[[1]])*abs(est[[1]])^(transf)
        est[[2]] <- sign(est[[2]])*abs(est[[2]])^(transf)
    }
    
    thin1.ind <- seq(1, length(ev[[1]]), by=thin)
    thin2.ind <- seq(1, length(ev[[2]]), by=thin)
    evx <- ev[[1]][thin1.ind]
    evy <- ev[[2]][thin2.ind]
    fx <- est[[1]][thin1.ind, thin2.ind]
    fy <- est[[2]][thin1.ind, thin2.ind]
    if (neg.grad) { fx <- -fx; fy <- -fy }
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
   
    if (!add) plot(fhat, abs.cont=c(0,0), col="transparent", xlab=xlab, ylab=ylab)
    grid.xy <- pracma::meshgrid(evy, evx)
    x0 <- grid.xy$Y
    y0 <- grid.xy$X
    
    ## default scale factor - arrows should exceed a bin
    arrow.len <- sqrt(fx^2 + fy^2) 
    bin.diag <- min(diff(evx), diff(evy)) 
    if (missing(scale)) scale <- bin.diag/max(arrow.len)
    
    ## remove `zero-length' arrows i.e. length < 1e-3 inches
    ## adapted from https://stackoverflow.com/questions/52689959/how-small-is-a-zero-length-arrow/52690054
    units <- par(c('usr', 'pin'))
    x2in <- with(units, pin[1L]/diff(usr[1:2]))
    y2in <- with(units, pin[2L]/diff(usr[3:4]))
    arrow2in <- sqrt((x2in*scale*fx)^2 + (y2in*scale*fy)^2)
    
    arrowint <- sort(c(min(arrow2in)-0.1*abs(max(arrow2in)), seq(1e-3, length, length=5), max(arrow2in)+0.1*abs(max(arrow2in))))
    nsmall.f <- cut(arrow2in, arrowint, labels=FALSE)
    nsf.lab <- unique(nsmall.f[nsmall.f>1])
    
    for (i in nsf.lab)
    {
        nsf <- which(nsmall.f==i)
        if (length(nsf)>0) pracma::quiver(x=x0[nsf], y=y0[nsf], u=fx[nsf], v=fy[nsf], col=col, scale=scale, length=arrowint[i], ...)
    }
}

## contourLevels method 

contourLevels.kdde <- function(x, prob, cont, nlevels=5, approx=TRUE, which.deriv.ind=1, ...)
{ 
    fhat <- x
    if (is.vector(fhat$x))
    {
        d <- 1; n <- length(fhat$x)
        if (!is.null(fhat$deriv.order))  
        {
          fhat.temp <- fhat 
          fhat.temp$deriv.ind <-fhat.temp$deriv.ind[which.deriv.ind]
          fhat <- fhat.temp
        }    
    }
    else
    {
        d <- ncol(fhat$x); n <-nrow(fhat$x)
        if (!is.matrix(fhat$x)) fhat$x <- as.matrix(fhat$x)

        if (!is.null(fhat$deriv.order))  
        {
            fhat.temp <- fhat 
            fhat.temp$estimate <- fhat.temp$estimate[[which.deriv.ind]]
            fhat.temp$deriv.ind <-fhat.temp$deriv.ind[which.deriv.ind,]
            fhat <- fhat.temp
        }
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
            dobs <- predict.kde(fhat, x=fhat$x)
        else
            dobs <- kdde(x=fhat$x, H=fhat$H, eval.points=fhat$x, w=w, deriv.order=fhat$deriv.order)$estimate[,which.deriv.ind] 
            
        if (is.null(fhat$deriv.order))
        {
            if (!missing(prob) & missing(cont)) hts <- quantile(dobs[dobs>=0], prob=prob)
            if (missing(prob) & !missing(cont)) hts <- quantile(dobs[dobs>=0], prob=(100-cont)/100)
        }
        else
        {
            if (!missing(prob) & missing(cont)) hts <- rbind(-quantile(abs(dobs[dobs<0]), prob=prob), quantile(dobs[dobs>=0], prob=prob))
            if (missing(prob) & !missing(cont)) hts <- rbind(-quantile(abs(dobs[dobs<0]), prob=(100-cont)/100), quantile(dobs[dobs>=0], prob=(100-cont)/100))
        }
    }
  
    return(hts)
}

## predict method for KDDE 

predict.kdde <- function(object, ..., x)
{
    fhat <- object
    if (is.vector(fhat$H)) d <- 1 else d <- ncol(fhat$H)
    if (d==1) n <- length(x)
    else
    {
        if (is.vector(x)) x <- matrix(x, nrow=1)
        else x <- as.matrix(x)
        n <- nrow(x)
    }

    if (!is.null(fhat$deriv.ind))
    {
        if (is.vector(fhat$deriv.ind)) pk.mat <- predict.kde(fhat, x=x, ...)
        else
        {
            nd <- nrow(fhat$deriv.ind)
            pk.mat <- matrix(0, ncol=nd, nrow=n)
            for (i in 1:nd)
            {
                fhat.temp <- fhat
                fhat.temp$estimate <- fhat$estimate[[i]]
                pk.mat[,i] <- predict.kde(fhat.temp, x=x, ...)
            }
        }
    }
    else
        pk.mat <- predict.kde(fhat, x=x, ...)
    
    return(drop(pk.mat))
}

######################################################################
## Summary kernel curvature 
######################################################################

kcurv <- function(fhat, compute.cont=TRUE)
{
    fhat.curv <- fhat
    if (is.vector(fhat$H)) d <- 1 else d <- ncol(fhat$H)
    if (fhat$deriv.order!=2) stop("Requires output from kdde(, deriv.order=2).")

    if (d==1)
    {
        Hessian.det <- fhat$estimate 
        local.mode <- fhat$estimate <0
        fhat.curv$estimate <- local.mode*abs(Hessian.det)
    }
    else if (d>1)
    {
        fhat.est <- sapply(fhat$estimate, as.vector)
        Hessian.det <- sapply(seq(1,nrow(fhat.est)), function(i) { det(invvec(fhat.est[i,])) })
        Hessian.eigen <- lapply(lapply(seq(1,nrow(fhat.est)), function(i) { invvec(fhat.est[i,]) }), eigen, only.values=TRUE)
        Hessian.eigen <- t(sapply(Hessian.eigen, getElement, "values"))
        local.mode <- apply(Hessian.eigen <= 0, 1, all)
        fhat.curv$estimate <- local.mode*array(abs(Hessian.det), dim=dim(fhat$estimate[[1]]))
    }
    
    fhat.curv$deriv.order <- NULL
    fhat.curv$deriv.ind <- NULL
    if (compute.cont)
    {
        fhat.temp <- fhat.curv
        fhat.temp$x <- fhat.curv$x[predict(fhat.curv, x=fhat.curv$x)>0,]
        fhat.temp$estimate <- fhat.temp$estimate
        fhat.curv$cont <- contourLevels(fhat.temp, cont=1:99)
    }
    class(fhat.curv) <- "kde"
    fhat.curv$type <- "kcurv"
    
    return(fhat.curv)
}
