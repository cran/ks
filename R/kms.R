###############################################################################
## Kernel mean shift
###############################################################################

kms <- function(x, y, H, max.iter=400, tol.iter, tol.clust, min.clust.size, merge=TRUE, keep.path=FALSE, verbose=FALSE)
{
    n <- nrow(x)
    d <- ncol(x)
    if (missing(tol.iter)) tol.iter <- 1e-3*min(apply(x, 2, IQR)) ## mean(apply(apply(x, 2, range), 2, diff))
    if (missing(tol.clust)) tol.clust <- 1e-2*max(apply(x, 2, IQR)) ## mean(apply(apply(x, 2, range), 2, diff))  
    if (missing(y)) y <- x
    if (missing(min.clust.size)) min.clust.size <- round(1e-2*nrow(y),0) 
    if (missing(H)) H <- Hpi(x, deriv.order=1, binned=default.bflag(d=d, n=n), nstage=2-(d>2))
    Hinv <- chol2inv(chol(H))
    if (is.vector(y)) y <- matrix(y, nrow=1)
   
    ## mean shift iterations
    n.seq <- block.indices(n, nrow(y), d=d, r=0, diff=FALSE, block.limit=3e6)
    if (verbose) pb <- txtProgressBar() 
    ms <- list()
    i <- 1
    if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    ms <- kms.base(x=x, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, tol.clust=tol.clust, Hinv=Hinv, verbose=verbose, max.iter=max.iter)
    
    if (length(n.seq)>2)
    {
        for (i in 2:(length(n.seq)-1))
        {  
            if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
            ms.temp <- kms.base(x=x, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, tol.clust=tol.clust, Hinv=Hinv, verbose=verbose, max.iter=max.iter)

            ms$y <- rbind(ms$y, ms.temp$y)
            ms$end.points <- rbind(ms$end.points, ms.temp$end.points)
            ms$label <- c(ms$label, ms.temp$label + max(ms$label))
            ms$mode <- rbind(ms$mode, ms.temp$mode)
            ms$nclust <- ms$nclust + ms.temp$nclust
            ms$nclust.table <- table(ms$label)
            ms$path <- c(ms$path, ms.temp$path)

            ## merge clusters which are closer than tol.clust distance
            ms <- ms.merge.dist(ms=ms, tol=tol.clust, verbose=FALSE)
        }     
    }
    if (verbose) close(pb)
    path.temp <- ms$path
    ms$path <- NULL
    ms$tol.iter <- tol.iter
    ms$tol.clust <- tol.clust
    ms$min.clust.size <- min.clust.size
    ms$names <- parse.name(x)
    if (keep.path) ms$path <- path.temp
    
    ## merge clusters which are smaller than min.clust.size
    if (merge) ms <- ms.merge.num(ms, min.clust.size=min.clust.size, verbose=verbose)

    return(ms)
}

kms.base <- function(x, H, Hinv, y, max.iter, tol.iter, tol.clust, verbose=FALSE)
{
    ## mean shift iterations
    ## original implementation J.E. Chacon (2013)
    ## modifications T.D. (2014)
    
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

    disp.ind <- head(sample(1:nrow(y)), n=min(100,nrow(y)))
    while (eps > tol.iter & i< max.iter)
    {
        y.curr <- y.update
        yHinvy <- t(rowSums(y.curr%*%Hinv *y.curr))
        Mah <- apply(yHinvy, 2, "+", xHinvx) - 2*xHinv %*% t(y.curr)
        w <- exp(-Mah/2)
        denom <- colSums(w)
        num <- t(w)%*%x
        denom[denom<=1e-3*tol.iter] <- 1e-3*tol.iter
        mean.shift.H <- num/denom
        y.update <- mean.shift.H                             
        y.update.list <- split(y.update, row(y.update), drop=FALSE)
        y.path <- mapply(rbind, y.path, y.update.list, SIMPLIFY=FALSE)
        eps <- max(sqrt(rowSums((y.curr-y.update)^2)))
    
        if (verbose>1)
        {
            y.range <- apply(y, 2, range)
            if (d==2) plot(y.update[disp.ind,], col=1, xlim=y.range[,1], ylim=y.range[,2], xlab="x", ylab="y")
            else pairs(y.update[disp.ind,], col=1)
        }
        i <- i+1
    }
    ms.endpt <- t(sapply(y.path, tail, n=1, SIMPLIFY=FALSE))

    ## extract cluster centres
    mode.tree <- hclust(dist(ms.endpt))    
    clust.label <- cutree(mode.tree, h=tol.clust)
    nclust <- length(unique(clust.label))
    mode.val <- by(ms.endpt, INDICES=clust.label, FUN=colMeans)
    mode.val  <- t(sapply(mode.val, FUN=identity))
    colnames(mode.val) <- colnames(x)
    rownames(mode.val) <- NULL
    nclust.table <- table(clust.label, dnn="")
    
    ms <- list(x=x, y=y, end.points=ms.endpt, H=H, label=clust.label, nclust=nclust, nclust.table=nclust.table, mode=mode.val, path=y.path)
    class(ms) <- "kms"

    return(ms)
}

## merge classes in 'label' into a single class
## label is a list of vector of class labels
## i.e. all classes in label[[j]] merged into new class j

ms.merge.label <- function(ms, label, verbose=FALSE)
{
    ms.merge <- ms
    for (i in 1:length(label))
    {
        labeli <- label[[i]]
       
        if (length(labeli)>1)
        {    
            merge.label <- min(labeli)
            ms.merge$label[ms$label %in% labeli] <- merge.label
            mode.label <- ms$mode[labeli[round(length(labeli)/2,0)],]
            for (j in labeli) ms.merge$mode[j,] <- mode.label
        }
    }

    ms.merge$mode <- unique(ms.merge$mode)
    ms.merge$nclust <- nrow(ms.merge$mode)
    ms.merge$label <- as.factor(ms.merge$label)
    levels(ms.merge$label) <- 1:ms.merge$nclust
    ms.merge$label <- as.numeric(ms.merge$label)
    ms.merge$nclust.table <- table(ms.merge$label)

    if (verbose) cat("Current clusters:", ms.merge$nclust.table, "\n")
    
    return(ms.merge)  
}
   
## merge mean shift clusters based on distance threshold

ms.merge.dist <- function(ms, tol, verbose)
{   
    if (missing(tol)) tol <- 1e-1*min(apply(ms$x, 2, IQR))
    mode.tree <- hclust(dist(ms$mode))    
    merge.label <- cutree(mode.tree, h=tol)
    merge.label <- split(1:ms$nclust, merge.label, drop=FALSE)

    ## create list where each element is a vector of cluster labels
    ## to be merged into a single cluster
    ms.temp <- ms.merge.label(ms=ms, label=merge.label, verbose=verbose)
    
    return(ms.temp)
}

## merge mean shift clusters based on min cluster size

ms.merge.num <- function(ms, min.clust.size, verbose=FALSE)
{
    if (missing(min.clust.size)) min.clust.size <- round(1e-2*nrow(ms$y),0)
    min.clust.size <- round(min.clust.size, 0)

    if (any(ms$nclust.table<=min.clust.size))
    {
        if (verbose) cat("Min cluster size merging begins. Min size = ", min.clust.size, "\n")

        ms.temp <- ms
        while(any(ms.temp$nclust.table<=min.clust.size) & ms.temp$nclust>1)
        {
            nclust.table <- table(ms.temp$label)
            small.clust.ind <- which.min(nclust.table)
            
            if (nclust.table[small.clust.ind] <= min.clust.size)
            {
                nearest.clust.ind <- FNN::get.knnx(ms.temp$mode, ms.temp$mode, k=2)$nn.index[small.clust.ind,2]

                merge.label <- ms$label
                merge.label[merge.label==small.clust.ind] <- nearest.clust.ind 
                merge.label <- split(ms$label, merge.label, drop=FALSE)
                merge.label <- lapply(merge.label, unique)
                ms.temp <- ms.merge.label(ms=ms.temp, label=merge.label, verbose=verbose)
            }
        }
        ms <- ms.temp
        ms$min.clust.size <- min.clust.size
        
        if (verbose) cat("Min cluster size merging ends.\n\n")
    }
    if (verbose) { cat("Final clusters:\n"); summary(ms) }

    return(ms)
}

######################################################################
## Cluster partition for 2D kernel mean shift
#####################################################################

kms.part <- function(x, H, xmin, xmax, gridsize, verbose=FALSE, ...)
{
    if (missing(H)) H <- Hpi(x, deriv.order=1, binned=TRUE)
    tol <- 5
    tol.H <-  tol * diag(H)
    if (missing(xmin)) xmin <- apply(x, 2, min) - tol.H
    if (missing(xmax)) xmax <- apply(x, 2, max) + tol.H
    if (missing(gridsize)) gridsize <- default.gridsize(2) 
    xx <- seq(xmin[1], xmax[1], length = gridsize[1])
    yy <- seq(xmin[2], xmax[2], length = gridsize[2])
    xy <- expand.grid(xx, yy)

    xy.kms <- kms(x=x, y=xy, H=H, verbose=verbose, ...)
    xy.lab <- array(xy.kms$label, dim=gridsize)
    
    fhat <- kde(x=x, binned=TRUE, xmin=xmin, xmax=xmax, bgridsize=gridsize)
    fhat$estimate <- xy.lab

    fhat <- c(fhat, xy.kms[c("end.points", "label", "mode", "nclust", "nclust.table", "min.clust.size", "tol.iter", "tol.clust")])
    class(fhat) <- "kde.part"
    
    return(fhat)
}

plot.kde.part <- function(x, display="filled.contour", col, col.fun, alpha=1, add=FALSE, ...)
{
    clev <- sort(unique(as.vector(x$estimate)))
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Dark2", alpha=alpha) }
    if (missing(col)) col <- col.fun(length(clev))
    
    for (i in 1:length(clev))
    {
        xtemp <- x
        xtemp$estimate <- x$estimate==clev[i]
        plot.kde(xtemp, display=display, col=c("transparent", col[i]), add=add | i>1, abs.cont=0.5, drawlabels=FALSE, alpha=alpha, ...)
    }
}

#############################################################################
## S3 methods for KMS objects
#############################################################################

## summary method

summary.kms <- function(object, ...)
{
    cat("Number of clusters =", object$nclust, "\n")
    cat("Cluster label table =", object$nclust.table, "\n")
    cat("Cluster modes =\n")
    print(as.data.frame(object$mode), ...)
}

## plot method

plot.kms <- function(x, display="splom", col, col.fun, alpha=1, xlab, ylab, zlab, theta=-30, phi=40, add=FALSE, ...)
{
    disp <- match.arg(display, c("splom", "plot3D", "rgl"))
    if (is.vector(x$H)) d <- 1 else d <- ncol(x$H)
    if (missing(col.fun)) col.fun <- function(n) { hcl.colors(n, palette="Set2", alpha=alpha) }
    if (missing(col)) col <- col.fun(length(unique(x$label)))
    col <- transparency.col(col, alpha=alpha)
    if (d==1) stop("kms plot not yet implemented")
    else if (d==2)
    {
        if (missing(xlab)) xlab <- x$names[1]
        if (missing(ylab)) ylab <- x$names[2]
        if (!add) plot(x$x, col=col[x$label], xlab=xlab, ylab=ylab, ...)
        else points(x$x, col=col[x$label], ...)
    }
    else if (d==3 & disp %in% c("plot3D", "rgl"))
    {
        if (missing(xlab)) xlab <- x$names[1]
        if (missing(ylab)) ylab <- x$names[2]
        if (missing(zlab)) zlab <- x$names[3]
    
        if (disp=="plot3D")
        {
            if (!requireNamespace("plot3D", quietly=TRUE)) stop("Install the plot3D package as it is required.", call.=FALSE)
            if (!add) plot3D::points3D(x$x[,1], x$x[,2], x$x[,3], col=1, cex=0, add=add, theta=theta, phi=phi, d=4, colkey=FALSE, xlab=xlab, ylab=ylab, zlab=zlab, ticktype="detailed", bty="f", ...)
        
            for (i in 1:length(col))
                plot3D::points3D(x$x[x$label==i,1], x$x[x$label==i,2], x$x[x$label==i,3], col=col[i], add=TRUE, ...)
        }
        else if (disp=="rgl")
        {
            ## suggestions from Viktor Petukhov 08/03/2018
            if (!requireNamespace("rgl", quietly=TRUE)) stop("Install the rgl package as it is required.", call.=FALSE)
            if (!add) rgl::plot3d(x$x, col=col[x$label], alpha=alpha, ...)
            else rgl::points3d(x$x, col=col[x$label], alpha=alpha, ...) 
        }
    }
    else if (d>=3)
    {
        pairs(x$x, col=col[x$label], ...)
    }    
}
