######################################################################
## Kernel density ridge estimation for 2D/3D data
#####################################################################

kdr <- function(x, y, H, p=1, max.iter=400, tol.iter, segment=TRUE, k, kmax, min.seg.size, keep.path=FALSE, gridsize, xmin, xmax, binned, bgridsize, w, fhat, density.cutoff, pre=TRUE, verbose=FALSE)
{
    ## default values 
    xnames <- parse.name(x)
    x <- as.matrix(x)
    x.orig <- x
    if (pre) 
    { 
        S12 <- diag(apply(x.orig, 2, sd))  
        Sinv12 <- matrix.pow(S12,-1)
        x <- pre.scale(x)
        if (!missing(xmin)) xmin <- xmin %*% Sinv12
        if (!missing(xmax)) xmax <- xmax %*% Sinv12
        rescale <- function(x) { as.matrix(x) %*% S12 }   
    }
    
    ksd <- ks.defaults(x=x, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize

    ## default bandwidth
    if (missing(H)) H <- Hpi(x=x, nstage=2-(d>2), binned=binned, deriv.order=2, verbose=verbose)
    Hinv <- chol2inv(chol(H))
    tol <- 3.7
    tol.H <-  tol * diag(H)

    if (missing(xmin)) xmin <- apply(x, 2, min) - tol.H 
    if (missing(xmax)) xmax <- apply(x, 2, max) + tol.H
    if (missing(tol.iter)) tol.iter <- 1e-3*min(apply(x, 2, IQR))
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
    else { y <- as.matrix(y); if (pre) y <- y %*% Sinv12 }
   
    if (is.vector(y)) y <- matrix(y, nrow=1)
    if (missing(min.seg.size)) min.seg.size <- round(1e-3*nrow(y), 0)
    
    ## exclude low density regions from ridge search 
    if (missing(fhat)) fhat <- kde(x=x, w=w, binned=binned)
    if (missing(density.cutoff)) density.cutoff <- contourLevels(fhat, cont=99)
    y.ind <- predict(fhat, x=y)>density.cutoff
    y <- y[y.ind,]
   
    fhat2 <- kdde(x=x, H=H, deriv.order=2, xmin=xmin, xmax=xmax, binned=binned, bgridsize=bgridsize, gridsize=gridsize, w=w, verbose=verbose) 
    
    ## projected gradient mean shift iterations
    n.seq <- block.indices(n, nrow(y), d=d, r=0, diff=FALSE)#, block.limit=1e6)
    if (verbose) pb <- txtProgressBar() 
    pc <- list()
    
    i <- 1
    if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
    pc <- kdr.base(x=x, fhat2=fhat2, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, Hinv=Hinv, verbose=verbose, max.iter=max.iter, p=p)
    
    if (pre)
    {
        pc[c("x","y","end.points")] <- lapply(pc[c("x","y","end.points")], rescale)
        pc[["path"]] <- lapply(pc[["path"]], rescale)
    }
    
    if (length(n.seq)>2)
    {
        for (i in 2:(length(n.seq)-1))
        {  
            if (verbose) setTxtProgressBar(pb, i/(length(n.seq)-1))
            pc.temp <- kdr.base(x=x, fhat2=fhat2, y=y[n.seq[i]:(n.seq[i+1]-1),], H=H, tol.iter=tol.iter, Hinv=Hinv, verbose=verbose, max.iter=max.iter, p=p)

            if (pre)
            {
                pc.temp[c("y","end.points")] <- lapply(pc.temp[c("y","end.points")], rescale)
                pc.temp[["path"]] <- lapply(pc.temp[["path"]], rescale)
            }
            pc$y <- rbind(pc$y, pc.temp$y)
            pc$end.points <- rbind(pc$end.points, pc.temp$end.points)
            pc$path <- c(pc$path, pc.temp$path)
        }
    }
    if (verbose) close(pb)

    ## remove short segments for p=1
    if (p==1)
    {     
        tol.seg <- 1e-2*max(apply(x.orig, 2, IQR))
        pc.dendo <- hclust(dist(pc$end.points), method="single")
        pc.label <- cutree(pc.dendo, h=tol.seg)
        pc.label.ind <- pc.label %in% which(table(pc.label)>min.seg.size)
        pc$y <- pc$y[pc.label.ind,]
        pc$end.points <- pc$end.points[pc.label.ind,]
        pc$path <- pc$path[pc.label.ind]
    }
    pc$H <- H
    pc$names <- xnames
    if (pre) pc$H <- S12 %*% pc$H %*% S12
    
    if (segment) pc <- kdr.segment(x=pc, k=k, kmax=kmax, min.seg.size=min.seg.size, verbose=verbose)
    else pc$end.points <- data.frame(pc$end.points, segment=1L)
     
    ## put paths as last element in list
    path.temp <- pc$path
    pc$path <- NULL
    pc$tol.iter <- tol.iter
    pc$min.seg.size <- min.seg.size
    pc$binned <- binned
    pc$names <- xnames
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
    
    pc <- list(x=x, y=y, end.points=pc.endpt, path=y.path, type="kdr")
    class(pc) <- "kdr"
    
    return(pc)
}

## create segment of KDR filaments
## x = output from kdr

kdr.segment <- function(x, k, kmax, min.seg.size, verbose=FALSE)
{
    ep <- x$end.points
    if (any(names(ep) %in% "segment")) ep <- x$end.points[,-which(names(ep)=="segment")]
    if (missing(min.seg.size)) min.seg.size <- x$min.seg.size
    hc <- hclust(dist(ep), method="single")
    
    if (missing(kmax)) kmax <- 30
    kmax <- min(kmax, nrow(ep))
    if (missing(k))
    {
        if (verbose) pb <- txtProgressBar() 
         
        clust.ind <- rep(0, kmax)
        for (i in 1:kmax)
        {
            if (verbose) setTxtProgressBar(pb, i/kmax)
            clust.ind[i] <- clust.crit(hc=hc, x=as.matrix(ep), k=i, min.seg.size=min.seg.size)
        }
        if (verbose) close(pb)
        clust.ind[is.na(clust.ind)] <- 0
        kopt <- which.max(clust.ind) 
    }
    else kopt <- k
   
    ep <- data.frame(ep, segment=as.integer(cutree(hc, k=kopt)))
    label <- ep$segment
    tlabel <- as.integer(names(table(label))[table(label) > min.seg.size])
    ep <- ep[label %in% tlabel,]
    ep$segment <- factor(ep$segment, labels=1:length(unique(ep$segment)))
    ep$segment <- as.integer(levels(ep$segment))[ep$segment]
    
    ## re-order KDR segments into 'reasonable' linestring order
    ## experimental 
    j <- 1
    for (i in unique(ep$segment))
    {
        ep.temp <- as.matrix(ep[ep$segment==i,-ncol(ep)])
        ep.temp <- data.frame(chain.knnx(ep.temp, k1=1, k2=1), segment=i) 
        if (j==1) ep.ord <- ep.temp else ep.ord <- rbind(ep.ord, ep.temp) 
        j <- j+1
    }
    names(ep.ord) <- c(x$names, "segment")
    rownames(ep.ord) <- NULL
    x$end.points <- ep.ord
    
    x$min.seg.size <- min.seg.size
    x$k <- kopt
    if (exists("clust.ind")) x$clust.ind <- clust.ind
    
    return(x)
}

## rbind nearest neighbour of y from x to y

add.knnx <- function(x, y, k=1)
{
    xynn <- FNN::get.knnx(x, y, k=k)
    y <- rbind(y, x[xynn$nn.index,])
    d <- ncol(x)
    xy <- list(x=matrix(x[-xynn$nn.index,],ncol=d), y=y)
    
    return(xy)      
}

## arrange points in KDR to form a "reasonable" linestring 

chain.knnx <- function(x, k1=1, k2=5)
{
    ## concatenate the nearest neighbours in a chain
    ## start with first point in x
    if (!is.matrix(x)) x <- as.matrix(x)
    d <- ncol(x)
    if (nrow(x)>1)
    {
        x.ord.list <- add.knnx(x=matrix(x[-1,], ncol=d), y=matrix(x[1,], ncol=d), k=k1)
        x.ord <- x.ord.list$y

        while (nrow(x.ord.list$x)>0)
        {
            y.temp <- matrix(apply(as.matrix(tail(x.ord.list$y, n=k2)), 2, mean), ncol=d)
            colnames(y.temp) <- names(x.ord.list$x)
            x.ord.list.temp <- add.knnx(x=x.ord.list$x, y=y.temp, k=k1)  
            x.ord <- rbind(x.ord, matrix(tail(x.ord.list.temp$y,n=1),ncol=d))
            x.ord.list <- x.ord.list.temp
        }
        
        ## decide which permutation is "best" linestring
        ## break at max discontinuity
        ind <- which.max(rowSums((head(x.ord, n=-1)-tail(x.ord,n=-1))^2))
        ind1 <- c(1:ind, (ind+1):nrow(x.ord))
        ind2 <- c(1:ind, rev((ind+1):nrow(x.ord)))
        ind3 <- c(rev(1:ind), (ind+1):nrow(x.ord))
        ind4 <- c(rev(1:ind), rev((ind+1):nrow(x.ord)))
        x.ord1 <- x.ord[ind1,] 
        x.ord2 <- x.ord[ind2,]
        x.ord3 <- x.ord[ind3,] 
        x.ord4 <- x.ord[ind4,]
        x.ord.dist <- rep(0,4)
        x.ord.dist[1] <- sum(rowSums((head(x.ord1, n=-1)-tail(x.ord1,n=-1))^2))
        x.ord.dist[2] <- sum(rowSums((head(x.ord2, n=-1)-tail(x.ord2,n=-1))^2))
        x.ord.dist[3] <- sum(rowSums((head(x.ord3, n=-1)-tail(x.ord3,n=-1))^2))
        x.ord.dist[4] <- sum(rowSums((head(x.ord4, n=-1)-tail(x.ord4,n=-1))^2))
        x.ord <- get(paste0("x.ord", which.min(x.ord.dist)))
    }
    else x.ord <- x
    
    return(x.ord)
}

## Calinski-Harabasz clustering criterion for hierarchical clustering object

clust.crit <- function(hc, x, k, min.seg.size=1)
{
    label <- cutree(hc, k=k)
    tlabel <- as.integer(names(table(label))[table(label) > min.seg.size])
    tlabel.ind <- label %in% tlabel
    cc <- fpc.calinhara(x=x, clustering=label)
    
    return(cc)
}

## copied from fpc::calinhara 2020-09-18

fpc.calinhara <- function(x, clustering, cn = max(clustering)) 
{
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    cln <- rep(0, cn)
    W <- matrix(0, p, p)
    for (i in 1:cn) cln[i] <- sum(clustering == i)
    for (i in 1:cn) {
        clx <- x[clustering == i, ]
        cclx <- cov(as.matrix(clx))
        if (cln[i] < 2) 
            cclx <- 0
        W <- W + ((cln[i] - 1) * cclx)
    }
    S <- (n - 1) * cov(x)
    B <- S - W
    out <- (n - cn) * sum(diag(B))/((cn - 1) * sum(diag(W)))
    
    return(out)
}

#############################################################################
## S3 methods for KDR objects
#############################################################################

## plot method 

plot.kdr <- function(x, ...)
{ 
    fhat <- x
    d <- ncol(fhat$x)

    if (d==2) 
    {
        plotret <- plotkdr.2d(fhat, ...)
        invisible(plotret)
    }
    else if (d==3)
    {
        plotkdr.3d(fhat, ...)
        invisible()
    }
    else 
      stop ("Plot function only available for 2 or 3-d data")
}

plotkdr.2d <- function(x, add=FALSE, col, type="p", alpha=1, ...)
{
    xp <- x$end.points
    if (!any(names(xp) %in% "segment"))
    {
        if (missing(col)) col <- 6
        col <- transparency.col(col, alpha=alpha)
        if (!add) plot(xp, col=col, type=type, ...)
        else points(xp, col=col, ...)  
    }
    else
    {
        xps <- unique(xp$segment)
        if (missing(col)) col <- hcl.colors(n=length(xps), palette="Set2")
        if (length(col) < length(xps)) col <- rep(col, length(xps))
        col <- transparency.col(col, alpha=alpha)
        if (!add) plot(xp[,-ncol(xp)], col="transparent", ...)
        for (i in 1:length(xps)) lines(xp[xp$segment==xps[i],-ncol(xp)], col=col[i], type=type, ...)  
    }
}

plotkdr.3d <- function(x, display="plot3D", colors, col, col.fun, alphavec, size=3, cex=1, pch=1, theta=-30, phi=40, d=4, ticktype="detailed", add=FALSE, xlab, ylab, zlab, alpha=1, box=TRUE, axes=TRUE, type="p", ...)
{
    fhat <- x
    xp <- x$end.points
  
    if (missing(xlab)) xlab <- fhat$names[1]
    if (missing(ylab)) ylab <- fhat$names[2]
    if (missing(zlab)) zlab <- fhat$names[3]
    if (!any(names(xp) %in% "segment"))
    {
        if (missing(col)) col <- 6
    }
    else 
    {
        xps <- unique(xp$segment)
        if (missing(col)) col <- hcl.colors(n=length(xps), palette="Set2")
        if (length(col) < length(xps)) col <- rep(col, length(xps))
    }
    colors <- col
    disp <- match.arg(display, c("plot3D", "rgl")) 
	if (disp %in% "plot3D") 
	{
        if (!requireNamespace("plot3D", quietly=TRUE)) stop("Install the plot3D package as it is required.", call.=FALSE)
        if (!add) plot3D::scatter3D(x=xp[,1], y=xp[,2], z=xp[,3], add=add, theta=theta, phi=phi, d=d, type=type, xlab=xlab, ylab=ylab, zlab=zlab, ticktype=ticktype, type="n", col=NA, ...) 
        for (i in 1:length(xps))
            plot3D::scatter3D(x=xp[xp$segment==xps[i],1], y=xp[xp$segment==xps[i],2], z=xp[xp$segment==xps[i],3], cex=cex, col=col[i], add=TRUE, pch=pch, type=type, alpha=alpha, ...) 
    }
    else if (disp %in% "rgl")
    {
        if (!requireNamespace("rgl", quietly=TRUE)) stop("Install the rgl package as it is required.", call.=FALSE)
        for (i in 1:length(xps))
            rgl::plot3d(x=xp[xp$segment==xps[i],1], y=xp[xp$segment==xps[i],2], z=xp[xp$segment==xps[i],3], col=col[i], alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add | (i>1), box=box, axes=axes, type=type, size=size, ...)
    }
}
