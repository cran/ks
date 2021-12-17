###############################################################################
## Feature significance for ultivariate kernel density stimate 
###############################################################################

kfs <- function(x, H, h, deriv.order=2, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE, signif.level=0.05)
{
    r <- deriv.order

    ## default values 
    ksd <- ks.defaults(x=x, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n; w <- ksd$w
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize
       
    if (d==1)
    {
        if (missing(h)) h <- hpi(x=x, nstage=2, binned=default.bflag(d=d, n=n), deriv.order=r)

        ## KDDE for r=2
        fhatr <- kdde(x=x, h=h, deriv.order=r, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, deriv.vec=FALSE, verbose=verbose)

        ## KDE
        fhat <- kde(x=x, h=h, gridsize=gridsize, gridtype=gridtype, xmin=min(fhatr$eval.points), xmax=max(fhatr$eval.points), binned=binned, bgridsize=bgridsize, positive=positive, adj.positive=adj.positive, w=w)

     
        fhat.est <- as.vector(fhat$estimate)
        fhatr.est <- as.vector(fhatr$estimate)

        RDrK <- (-1)^r*psins.1d(r=2*r, sigma=1)
        fhatr.Sigma <- n^(-1)* h^(-2*r-1)*RDrK
        fhatr.Sigma12 <- sqrt(fhatr.Sigma)
        fhatr.est <- fhatr.est/fhatr.Sigma12

        local.mode <- fhatr.est <= 0
        fhatr.wald <- fhatr.est^2/abs(fhat.est)
        fhatr.wald[is.infinite(fhatr.wald)] <- max(fhatr.wald[!is.infinite(fhatr.wald)])
        gs <- length(fhat$estimate)        
    }
    else if (d>1)
    {
        if (missing(H)) H <- Hpi(x=x, nstage=2-(d>2), binned=default.bflag(d=d, n=n), deriv.order=r, verbose=verbose)
    
        ## KDDE for r=2
        fhatr <- kdde(x=x, H=H, deriv.order=r, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, deriv.vec=FALSE, verbose=verbose)

        ## KDE
        fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, verbose=verbose)
                
        fhat.est <- as.vector(fhat$estimate)
        fhatr.est <- sapply(fhatr$estimate, as.vector)
    
        ## convert from vec to vech because vec'ed derivative
        ## contains repeated columns so its variance isn't invertible 
        Hinv <- chol2inv(chol(H))
        Hinv12 <- matrix.sqrt(Hinv)
        RDrK <- (-1)^r*invvec(psins(r=2*r, Sigma=diag(d)))
        dupld <- dupl(d)$d
        dupld.MPinv <- chol2inv(chol(t(dupld)%*% dupld)) %*% t(dupld)
        fhatr.Sigma.const <- n^(-1)*det(H)^(-1/2)* Kpow(Hinv12,r) %*% RDrK %*% Kpow(Hinv12,2)
        fhatr.Sigma.const <- dupld.MPinv %*% fhatr.Sigma.const %*% t(dupld.MPinv)      
        fhatr.Sigma.const12inv <- chol2inv(chol(matrix.sqrt(fhatr.Sigma.const)))
        fhatr.est <- fhatr.est %*% fhatr.Sigma.const12inv  

        ## all eigenvalues < 0 => local mode 
        fhatr.eigen <- lapply(lapply(seq(1,nrow(fhatr.est)), function(i) {invvech(fhatr.est[i,])}), eigen, only.values=TRUE)
        fhatr.eigen <- t(sapply(fhatr.eigen, getElement, "values"))
        local.mode <- apply(fhatr.eigen <= 0, 1, all)
        fhatr.wald <- apply(fhatr.est^2, 1, sum)/fhat.est
        
        gs <- dim(fhat$estimate)
    }

    ## Hochberg adjustment for sequential tests
    pval.wald <- 1 - pchisq(fhatr.wald, d*(d+1)/2)
    pval.wald[fhat.est<=contourLevels(fhat, cont=99) | !local.mode] <- NA    
    pval.wald.ord.index <- order(pval.wald)
    pval.wald.ord <- pval.wald[pval.wald.ord.index]
    num.test <- sum(!is.na(pval.wald.ord))
    
    if (num.test>=1) num.test.seq <- c(1:num.test, rep(NA, prod(gs) - num.test))
    else num.test.seq <- rep(NA, prod(gs))

    reject.nonzero <- (pval.wald.ord <= signif.level/(num.test + 1 - num.test.seq))
    reject.nonzero.ind <- which(reject.nonzero)
    
    ## reject null hypotheses indicated in reject.nonzero.ind
    if (d==1)
    {
        signif.wald <- array(0L, dim=gs)
        signif.wald[pval.wald.ord.index[reject.nonzero.ind]] <- 1L
    }
    else
    {
        signif.wald <- array(0L, dim=gs)
        signif.wald.index <- expand.grid(lapply(gs, seq_len))
        signif.wald[as.matrix(signif.wald.index[pval.wald.ord.index[reject.nonzero.ind],])] <- 1L
    }
    
    ## ESS = effective sample size
    ## ess <- n*fhat$estimate*dmvnorm.mixt(x=rep(0,d), mu=rep(0,d), Sigma=H, props=1)
    ## signif.ess <- ess >= 5

    fhatr$estimate <- signif.wald 
    fhatr$type <- "kfs"
    class(fhatr) <- "kfs"
    
    return(fhatr)
}

#############################################################################
## S3 methodfor KFS objects
#############################################################################

## plot method
plot.kfs <- function(x, display="filled.contour", col=7, colors, abs.cont, alpha=1, alphavec=0.4, add=FALSE, ...)
{
    fhatr <- x
    fhatr$deriv.order <- NULL
    class(fhatr) <- "kde"

    if (is.vector(fhatr$H)) d <- 1 else d <- ncol(fhatr$H)
    if (d==1)
    {
        fhat <- kde(x=fhatr$x, xmin=min(fhatr$eval.points), xmax=max(fhatr$eval.points), gridsize=length(fhatr$eval.points))
        plot(fhat, col="grey", add=add, ...)
        gridsize <- length(fhatr$estimate)

        estimate.rle <- rle(as.vector(fhatr$estimate))
        estimate.rle.cumsum <- rbind(cumsum(estimate.rle$lengths), estimate.rle$lengths, estimate.rle$values)
        col <- transparency.col(col, alpha=alpha)
        for (i in which(estimate.rle$values==1))
        {
            seg.ind <- 1:estimate.rle.cumsum[2,i]   
            if (i>1)
                seg.ind <- seg.ind + estimate.rle.cumsum[1,i-1] 
            lines(fhat$eval.points[seg.ind], fhat$estimate[seg.ind], col=col, ...)
        }
    }
    else if (d==2)
    {
        if (missing(abs.cont)) abs.cont <- 0.5
        disp1 <- match.arg(display, c("slice", "persp", "image", "filled.contour", "filled.contour2"))
        if (disp1=="filled.contour2") disp1 <- "filled.contour"
        col <- c("transparent",col)
        
        plot(fhatr, abs.cont=abs.cont, drawlabels=FALSE, col=col, add=add, display=display, alpha=alpha, ...)
    }
    else if (d==3)
    {
        if (missing(abs.cont)) abs.cont <- 0.25
        e1 <- try(match.arg(display, c("plot3D", "rgl")), silent=TRUE)
        if (class(e1) %in% "try-error") display <- "plot3D"
        if (!missing(colors)) col <- colors     
        plot(fhatr, abs.cont=abs.cont, col=col, colors=colors, alphavec=alphavec, add=add, display=display, ...)
    }
    invisible()
}

## predict method

predict.kfs <- function(object, ..., x)
{
    fhat <- predict.kde(object=object, ..., x=x, zero.flag=FALSE)
    fhat <- as.integer(fhat>=0.5)
   
    return(fhat)
}
