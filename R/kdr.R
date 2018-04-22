######################################################################
## Kernel density ridge estimation for 2D/3D data
#####################################################################

kdr <- function(x, y, H, p=1, max.iter=400, tol.iter, tol.seg, min.seg.size, keep.path=FALSE, gridsize, xmin, xmax, binned, bgridsize, w, fhat, density.cutoff, verbose=FALSE)
{
    ## default values 
    ksd <- ks.defaults(x=x, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize

    if (missing(tol.iter)) tol.iter <- 1e-3*min(apply(x, 2, IQR))
    if (missing(tol.seg)) tol.seg <- 1e-2*max(apply(x, 2, IQR))
    if (missing(H)) H <- Hpi(x=x, nstage=2-(d>2), binned=default.bflag(d=d, n=n), deriv.order=2, verbose=verbose)
    Hinv <- chol2inv(chol(H))
    tol <- 3.7
    tol.H <-  tol * diag(H)
    if (missing(xmin)) xmin <- apply(x, 2, min) - tol.H
    if (missing(xmax)) xmax <- apply(x, 2, max) + tol.H
    if (missing(y))
    {    
        xx <- seq(xmin[1], xmax[1], length = gridsize[1])
        yy <- seq(xmin[2], xmax[2], length = gridsize[2])
       
        if (d==2) y <- expand.grid(xx, yy)
        else if (d==3)
        {
            zz <- seq(xmin[3], xmax[3], length = gridsize[3])
            y <- expand.grid(xx, yy, zz)
        }
    }
    if (is.vector(y)) y <- matrix(y, nrow=1)
    if (missing(min.seg.size)) min.seg.size <- round(1e-3*nrow(y), 0)
    
    ## exclude low density regions from ridge search 
    if (missing(fhat)) fhat <- kde(x=x, w=w, binned=binned)
    if (missing(density.cutoff)) density.cutoff <- contourLevels(fhat, cont=99)
    y.ind <- predict(fhat, x=y)>density.cutoff
    y <- y[y.ind,]
   
    fhat2 <- kdde(x=x, H=H, deriv.order=2, xmin=xmin, xmax=xmax, binned=binned, bgridsize=bgridsize, gridsize=gridsize, w=w, verbose=verbose) 
    ## projected gradient mean shift iterations
    n.seq <- block.indices(n, nrow(y), d=d, r=0, diff=FALSE, block.limit=3e6)
    if (verbose) pb <- txtProgressBar() 
    pc <- list()
    i <- 1
    if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    pc <- kdr.base(x=x, fhat2=fhat2, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, Hinv=Hinv, verbose=verbose, max.iter=max.iter, p=p)
    
    if (length(n.seq)>2)
    {
        for (i in 2:(length(n.seq)-1))
        {  
            if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
            pc.temp <- kdr.base(x=x, fhat2=fhat2, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, Hinv=Hinv, verbose=verbose, max.iter=max.iter, p=p)

            pc$y <- rbind(pc$y, pc.temp$y)
            pc$end.points <- rbind(pc$end.points, pc.temp$end.points)
            pc$path <- c(pc$path, pc.temp$path)
        }
    }
    if (verbose) close(pb)
    
    ## remove short segments for p=1
    if (p==1)
    {
        pc.dendo <- hclust(dist(pc$end.points), method="single")
        pc.label <- cutree(pc.dendo, h=tol.seg)
        pc.label.ind <- pc.label %in% which(table(pc.label)>min.seg.size)
        pc$y <- pc$y[pc.label.ind,]
        pc$end.points <- pc$end.points[pc.label.ind,]
        pc$H <- H
        pc$path <- pc$path[pc.label.ind]
        pc$label <- factor(pc.label[pc.label.ind])
        levels(pc$label) <- 1:length(levels(pc$label))
        pc$label <- as.numeric(pc$label)
    }
    
    ## put paths as last element in list
    path.temp <- pc$path
    pc$path <- NULL
    pc$tol.iter <- tol.iter
    pc$tol.seg <- tol.seg
    pc$min.seg.size <- min.seg.size
    pc$binned <- binned
    pc$names <- parse.name(x)
    pc$w <- w
    if (keep.path) pc$path <- path.temp
    
    return(pc)
}


kdr.base <-function(x, fhat2, H, y, max.iter, tol.iter, p=1, verbose=FALSE, Hinv, ...)
{
    if (!is.matrix(x)) x <- as.matrix(x)
    if (!is.matrix(y)) y <- as.matrix(y)
    if (missing(Hinv)) Hinv <- chol2inv(chol(H))
    nx <- nrow(x)
    ny <- nrow(y)
    d <- ncol(y)
    y.path <- split(y, row(y), drop=FALSE)
    names(y.path) <- NULL
 
    xHinv <- x %*% Hinv
    xHinvx <- rowSums(xHinv*x)
    y.update <- y
    i <- 1
    eps <- max(sqrt(rowSums(y.update^2)))
    disp.ind <- head(sample(1:nrow(y)), n=min(1000,nrow(y)))

    while (eps > tol.iter & i < max.iter)
    {
	y.curr <- y.update
        yHinvy <- t(rowSums(y.curr%*%Hinv *y.curr))
        Mah <- apply(yHinvy, 2, "+", xHinvx) - 2*xHinv %*% t(y.curr)
        w <- exp(-Mah/2)
        denom <- colSums(w)
        num <- t(w)%*%x
        mean.shift.H <- num/denom - y.curr

        fhat2.y.curr <- predict(fhat2, x=y.curr)
        for (j in 1:ny)
        {
            Hessian <- invvec(fhat2.y.curr[j,])
            Hessian.svd <- eigen(Hessian, symmetric=TRUE)
            Up <- Hessian.svd$vectors[,tail(1:d,n=d-p)]
            mean.shift.H[j,] <- drop(Up %*% t(Up) %*% mean.shift.H[j,])
        }
        y.update <- y.curr + mean.shift.H                             
        y.update.list <- split(y.update, row(y.update), drop=FALSE)
        y.path <- mapply(rbind, y.path, y.update.list, SIMPLIFY=FALSE)
        eps <- max(sqrt(rowSums((y.curr-y.update)^2)))

        if (verbose>1)
        {
            if (d==2) plot(y.update[disp.ind,], col=1, xlab="x", ylab="y")
            else pairs(y.update[disp.ind,], col=1)
        }
        i <- i+1
    }
    pc.endpt <- t(sapply(y.path, tail, n=1, SIMPLIFY=FALSE))

    pc <- list(x=x, y=y, end.points=pc.endpt, path=y.path)
    class(pc) <- "kdr"
    
    return(pc)
}

