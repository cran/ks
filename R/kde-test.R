##############################################################################
## Test statistic for multivariate 2-sample test
##############################################################################

kde.test <- function(x1, x2, H1, H2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, binned=FALSE, bgridsize, verbose=FALSE)
{
    ## default values 
    ksd <- ks.defaults(x=x1, binned=binned, bgridsize=bgridsize)
    d <- ksd$d
    binned <- ksd$binned
    bgridsize <- ksd$bgridsize
    gridsize <- ksd$gridsize

    if (is.vector(x1) & is.vector(x2))
        return (kde.test.1d(x1=x1, x2=x2, h1=h1, h2=h2, psi1=psi1, psi2=psi2, var.fhat1=var.fhat1, var.fhat2=var.fhat2, binned=binned, bgridsize=bgridsize, verbose=verbose))

    if (!is.matrix(x1)) x1 <- as.matrix(x1)
    if (!is.matrix(x2)) x2 <- as.matrix(x2)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    d <- ncol(x1)
    K0 <- drop(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=0))
    S1 <- var(x1)
    S2 <- var(x2)

    ## default bandwidths 
    if (missing(H1)) H1 <- Hpi.kfe(x1, deriv.order=0, binned=default.bflag(d=d,n=n1), bgridsize=bgridsize)
    if (missing(H2)) H2 <- Hpi.kfe(x2, deriv.order=0, binned=default.bflag(d=d,n=n2), bgridsize=bgridsize)
    
    ## kernel estimation for components of test statistic
    if (missing(psi1)) psi1 <- Qr(x=x1, y=x1, Sigma=H1, verbose=verbose) 
    if (missing(psi2)) psi2 <- Qr(x=x2, y=x2, Sigma=H2, verbose=verbose)

    if (missing(var.fhat1))
    {
        H1.r1 <- Hns(x=x1, deriv.order=1)
        fhat1.r1 <- predict(kdde(x=x1, H=H1.r1, deriv.order=1), x=apply(x1, 2, mean))
        var.fhat1 <- drop(fhat1.r1 %*% S1 %*% fhat1.r1)
    }
    psi12 <- Qr(x=x1, y=x2, Sigma=H1, verbose=verbose) 

    if (missing(var.fhat2))
    {
        H2.r1 <- Hns(x=x2, deriv.order=1)
        fhat2.r1 <- predict(kdde(x=x2, H=H2.r1, deriv.order=1), x=apply(x2, 2, mean))
        var.fhat2 <- drop(fhat2.r1 %*% S2 %*% fhat2.r1)
        }
    psi21 <- Qr(x=x2, y=x1, Sigma=H2, verbose=verbose) 

    ## test statistic + its parameters
    T.hat <- drop(psi1 + psi2 - (psi12 + psi21))
    muT.hat <- (n1^(-1)*det(H1)^(-1/2) + n2^(-1)*det(H2)^(-1/2))*K0
    varT.hat <- 3*(n1*var.fhat1 + n2*var.fhat2)/(n1+n2) *(1/n1+1/n2) 

    zstat <- (T.hat-muT.hat)/sqrt(varT.hat)
    pval <- 1-pnorm(zstat)
    if (length(pval==0)>0) pval[pval==0] <- pnorm(-abs(zstat[pval==0]))

    val <- list(Tstat=T.hat, zstat=zstat, pvalue=pval, mean=muT.hat, var=varT.hat, var.fhat1=var.fhat1, var.fhat2=var.fhat2, n1=n1, n2=n2, H1=H1, H2=H2, psi1=psi1, psi12=psi12, psi21=psi21, psi2=psi2)

    return(val)
}     

kde.test.1d <- function(x1, x2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, binned=FALSE, bgridsize, verbose=FALSE)
{
    n1 <- length(x1)
    n2 <- length(x2)
    d <- 1
    K0 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=0)

    s1 <- sd(x1)
    s2 <- sd(x2)

    ## kernel estimation for components of test statistic
    if (missing(h1)) h1 <- hpi.kfe(x1, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize)
    if (missing(h2)) h2 <- hpi.kfe(x2, nstage=2, deriv.order=0, binned=binned, bgridsize=bgridsize)

    if (missing(psi1)) psi1 <- Qr.1d(x=x1, y=x1, sigma=h1, verbose=verbose)
    if (missing(psi2)) psi2 <- Qr.1d(x=x2, y=x2, sigma=h2, verbose=verbose)

    if (missing(var.fhat1))
    {
        h1.r1 <- hns(x=x1, deriv.order=1) 
        fhat1.r1 <- predict(kdde(x=x1, h=h1.r1, deriv.order=1), x=mean(x1))
        var.fhat1 <- fhat1.r1^2*s1^2
    }
    psi12 <- Qr.1d(x=x1, sigma=h1, y=x2, verbose=verbose)

    if (missing(var.fhat2))
    {
        h2.r1 <- hns(x=x2, deriv.order=1) 
        fhat2.r1 <- predict(kdde(x=x2, h=h2.r1, deriv.order=1), x=mean(x2))
        var.fhat2 <- fhat2.r1^2*s2^2
    }
    psi21 <- Qr.1d(x=x2, sigma=h2, y=x1, verbose=verbose)

    ## test statistic + its parameters
    T.hat <- drop(psi1 + psi2 - (psi12 + psi21))
    muT.hat <- ((n1*h1)^(-1) + (n2*h2)^(-1))*K0
    varT.hat <- 3*(n1*var.fhat1 + n2*var.fhat2)/(n1+n2) *(1/n1+1/n2) 
    zstat <- (T.hat-muT.hat)/sqrt(varT.hat)
    pval <- 1-pnorm(zstat)
    if (length(pval==0)>0) pval[pval==0] <- pnorm(-abs(zstat[pval==0])) 

    val <- list(Tstat=T.hat, zstat=zstat, pvalue=pval, mean=muT.hat, var=varT.hat, var.fhat1=var.fhat1, var.fhat2=var.fhat2, n1=n1, n2=n2, h1=h1, h2=h2, psi1=psi1, psi12=psi12, psi21=psi21, psi2=psi2)

    return(val)
}

###############################################################################
## Local kde test
###############################################################################

## multivariate local test

kde.local.test <- function(x1, x2, H1, H2, h1, h2, fhat1, fhat2, gridsize, binned, bgridsize, verbose=FALSE, supp=3.7, mean.adj=FALSE, signif.level=0.05, min.ESS, xmin, xmax)
{
    ## default values 
    ksd <- ks.defaults(x=x1, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d
    if (missing(binned)) binned <- ksd$binned
    if (missing(bgridsize)) bgridsize <- ksd$bgridsize
    if (missing(gridsize)) gridsize <- ksd$gridsize

    ## clip data to xmin, xmax grid limits for binned estimation
    grid.clip <- binned    
    if (grid.clip) 
    {
        if (!missing(xmax)) xmax <- xmax[1:d]
        if (!missing(xmin)) xmin <- xmin[1:d]
        x1 <- truncate.grid(x=x1, xmin=xmin, xmax=xmax)
        x2 <- truncate.grid(x=x2, xmin=xmin, xmax=xmax)
    }

    if (is.vector(x1) & is.vector(x2)) { return(kde.local.test.1d(x1=x1, x2=x2, h1=h1, h2=h2, fhat1=fhat1, fhat2=fhat2, gridsize=gridsize, binned=binned, bgridsize=bgridsize, verbose=verbose, supp=supp, mean.adj=mean.adj, xmin=xmin, xmax=xmax)) }

    if (missing(H1) & !missing(x1)) H1 <- Hpi(x=x1, deriv.order=0, binned=default.bflag(d=ncol(x1), n=nrow(x1)), bgridsize=bgridsize, verbose=verbose)
    if (missing(H2) & !missing(x2)) H2 <- Hpi(x=x2, deriv.order=0, binned=default.bflag(d=ncol(x2), n=nrow(x2)), bgridsize=bgridsize, verbose=verbose)

    if (!missing(x1) & !missing(x2))
    {
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        d <- ncol(x1)
        RK <- (4*pi)^(-d/2)

        xrange <- apply(rbind(x1,x2), 2, range)
        tol.H <- sqrt(prod(diag(matrix.sqrt(H1) %*% matrix.sqrt(H2)))^(1/d))
        if (missing(xmin)) xmin <- xrange[1,] - supp*tol.H 
        if (missing(xmax)) xmax <- xrange[2,] + supp*tol.H 
    }
    else
    {
        n1 <- nrow(fhat1$x)
        n2 <- nrow(fhat2$x)
        d <- ncol(fhat1$x)
        RK <- (4*pi)^(-d/2)
        H1 <- fhat1$H
        H2 <- fhat2$H 
    }

    ## kernel estimation for components of test statistic
    if (missing(fhat1)) fhat1 <- kde(x=x1, H=H1, binned=binned, gridsize=gridsize, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, verbose=verbose)
    if (missing(fhat2)) fhat2 <- kde(x=x2, H=H2, binned=binned, gridsize=gridsize, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, verbose=verbose)

    HD2fhat <- 0
    if (mean.adj)
    { 
        D2fhat1 <- kdde(x=x1, H=H1, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, deriv.order=2, verbose=verbose)
        D2fhat2 <- kdde(x=x2, H=H2, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, deriv.order=2, verbose=verbose)
        D2fhat <- list()
        vH1 <- vec(H1)
        vH2 <- vec(H2)
        for (j in 1:length(D2fhat1$estimate)) D2fhat$estimate[[j]] <- vH1[j]*D2fhat1$estimate[[j]] - vH2[j]*D2fhat2$estimate[[j]]
        for (j in 1:length(D2fhat)) HD2fhat <- HD2fhat + D2fhat$estimate[[j]]
    }

    fhat.diff <- fhat1$estimate - fhat2$estimate - 1/2*HD2fhat  
    var.fhat.diff <- (n1^(-1)*det(H1)^(-1/2)*fhat1$estimate + n2^(-1)*det(H2)^(-1/2)*fhat2$estimate)*RK
    X2 <- fhat.diff^2/var.fhat.diff
    pvalue <- 1 - pchisq(X2, 1)
    pvalue[is.na(pvalue)] <- 0

    ## apply Hochberg multiple test adjustment
    gridsize <- sapply(fhat1$eval.points, length)
    fhat.diff.signif <- hochberg.mult.test(pvalue=pvalue, gridsize=gridsize, signif.level=signif.level)    

    fhat.diff.pos <- fhat1
    fhat.diff.neg <- fhat2
    fhat.diff.pos$estimate <- fhat.diff.signif*(fhat.diff>0)
    fhat.diff.neg$estimate <- fhat.diff.signif*(fhat.diff<0)

    if (!missing(min.ESS))
    {
        ESS1 <- n1*fhat1$estimate/dmvnorm(rep(0,d), rep(0,d), H1)
        ESS2 <- n2*fhat2$estimate/dmvnorm(rep(0,d), rep(0,d), H2)
        ESS <- pmin(ESS1, ESS2) >= min.ESS
        fhat.diff.pos$estimate <- fhat.diff.pos$estimate*ESS
        fhat.diff.neg$estimate <- fhat.diff.neg$estimate*ESS
    }
    result <- list(fhat1=fhat1, fhat2=fhat2, X2=X2, pvalue=pvalue, fhat.diff=fhat.diff, mean.fhat.diff=HD2fhat, var.fhat.diff=var.fhat.diff, fhat.diff.pos=fhat.diff.pos, fhat.diff.neg=fhat.diff.neg, n1=n1, n2=n2, H1=H1, H2=H2, names=parse.name(x1))
    class(result) <- "kde.loctest"

    return(result)
}

## 1-d local test

kde.local.test.1d <- function(x1, x2, h1, h2, fhat1, fhat2, gridsize=gridsize, binned=FALSE, bgridsize, verbose=FALSE, supp=3.7, mean.adj=FALSE, signif.level=0.05, min.ESS, xmin, xmax)
{
    if (missing(h1) & !missing(x1)) h1 <- hpi(x=x1, nstage=2, binned=binned, bgridsize=bgridsize)
    if (missing(h2) & !missing(x2)) h2 <- hpi(x=x2, nstage=2, binned=binned, bgridsize=bgridsize)

    if (!missing(x1) & !missing(x2))
    {
        n1 <- length(x1)
        n2 <- length(x2)
        d <- 1 
        RK <- (4*pi)^(-d/2)
        xrange <- range(c(x1,x2))
        if (missing(xmin)) xmin <- xrange[1] - supp*sqrt(h1*h2)
        if (missing(xmax)) xmax <- xrange[2] + supp*sqrt(h1*h2)
    }
    else
    {
        n1 <- length(fhat1$x)
        n2 <- length(fhat2$x)
        d <- 1
        RK <- (4*pi)^(-d/2)
        h1 <- fhat1$h
        h2 <- fhat2$h 
    }

    ## kernel estimation for components of test statistic
    if (missing(fhat1)) fhat1 <- kde(x=x1, h=h1, gridsize=gridsize, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax)
    if (missing(fhat2)) fhat2 <- kde(x=x2, h=h2, gridsize=gridsize, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax)

    h2D2fhat <- 0
    if (mean.adj)
    { 
        D2fhat1 <- kdde(x=x1, h=h1, gridsize=gridsize, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, deriv.order=2, verbose=verbose)
        D2fhat2 <- kdde(x=x2, h=h2, gridsize=gridsize, binned=binned, bgridsize=bgridsize, supp=supp, xmin=xmin, xmax=xmax, deriv.order=2, verbose=verbose)
        h2D2fhat <- (h1^2*D2fhat1$estimate - h2^2*D2fhat2$estimate)
    }

    fhat.diff <- fhat1$estimate - fhat2$estimate - 1/2*h2D2fhat  
    var.fhat.diff <- ((n1*h1)^(-1)*fhat1$estimate + (n2*h2)^(-1)*fhat2$estimate)*RK

    X2 <- fhat.diff^2/var.fhat.diff
    pvalue <- 1 - pchisq(X2, 1)

    gridsize <- length(fhat1$eval.points)
    fhat.diff.signif <- hochberg.mult.test(pvalue=pvalue, gridsize=gridsize, signif.level=signif.level)    

    fhat.diff.pos <- fhat1
    fhat.diff.neg <- fhat2
    fhat.diff.pos$estimate <- fhat.diff.signif*(fhat.diff>0)
    fhat.diff.neg$estimate <- fhat.diff.signif*(fhat.diff<0)

    if (!missing(min.ESS))
    {
        ESS1 <- n1*fhat1$estimate/dnorm(0, 0, h1)
        ESS2 <- n2*fhat2$estimate/dnorm(0, 0, h2)
        ESS <- pmin(ESS1, ESS2) >= min.ESS
        fhat.diff.pos$estimate <- fhat.diff.pos$estimate*ESS
        fhat.diff.neg$estimate <- fhat.diff.neg$estimate*ESS
    }

    result <- list(fhat1=fhat1, fhat2=fhat2, chisq=X2, pvalue=pvalue, fhat.diff=fhat.diff, mean.fhat.diff=h2D2fhat, var.fhat.diff=var.fhat.diff, n1=n1, n2=n2, h1=h1, h2=h2, H1=h1^2, H2=h2^2, fhat.diff=fhat.diff, fhat.diff.pos=fhat.diff.pos, fhat.diff.neg=fhat.diff.neg, names=parse.name(x1))
    class(result) <- "kde.loctest"

    return(result)
}

## Hochberg (1988) adjustment for multiple correlated tests

hochberg.mult.test <- function(pvalue, gridsize, signif.level)
{
    pvalue.ord <- pvalue[order(pvalue)]
    num.test <- sum(!is.na(pvalue.ord))

    if (num.test>=1)
        num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else
        num.test.seq <- rep(NA, prod(gridsize))

    reject.nonzero <- ((pvalue.ord <= signif.level/(num.test + 1 - num.test.seq)) & (pvalue.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)

    ## p-value == 0 => reject null hypotheses automatically
    fhat.diff.signif <- array(FALSE, dim=gridsize)
    fhat.diff.signif[which(pvalue==0, arr.ind=TRUE)] <- TRUE

    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind 
    for (i in reject.nonzero.ind)
        fhat.diff.signif[which(pvalue==pvalue.ord[i], arr.ind=TRUE)] <- TRUE 

    return(fhat.diff.signif)    
}

#############################################################################
## S3 method for kde.loctest objects
#############################################################################

## plot methods

plot.kde.loctest <- function(x, ...)
{
    fhat <- x

    if (is.vector(fhat$fhat1$x)) plotkde.loctest.1d(fhat, ...)
    else
    {
        d <- ncol(fhat$fhat1$x)

        if (d==2) plotkde.loctest.2d(x, ...)
        else if (d==3) plotkde.loctest.3d(x, ...)
        else stop("Plot function only available for 1, 2 or 3-d data")   
    }
}

plotkde.loctest.1d <- function(x, lcol, col, add=FALSE, xlab, ylab, rugsize, add.legend=TRUE, pos.legend="topright", alpha=1, ...)
{
    if (missing(xlab)) xlab <- x$fhat.diff.pos$names[1]
    if (missing(ylab)) ylab <- expression("Density difference  "*f[1]-f[2])
    if (missing(col)) col <- hcl.colors(palette="Dark2",2)  
    col <- transparency.col(col, alpha=alpha)
    if (missing(lcol)) lcol <- 1
    if (!add) plot(x$fhat1$eval.points, x$fhat.diff, type="l", ylab=ylab, xlab=xlab, col=lcol, ...)
    else lines(x$fhat1$eval.points, x$fhat.diff, col=lcol, ...)

    plot.lim <- par()$usr
    if (missing(rugsize)) rugsize <- abs(plot.lim[4]-plot.lim[3])/50

    try(image(x$fhat.diff.pos$eval, c(plot.lim[3], plot.lim[3]+rugsize), cbind(x$fhat.diff.pos$estimate==1, x$fhat.diff.pos$estimate==1), level=0.5, add=TRUE, col=c("transparent", col[1]), ...))
    try(image(x$fhat.diff.neg$eval, c(plot.lim[3], plot.lim[3]+rugsize), cbind(x$fhat.diff.neg$estimate==1, x$fhat.diff.neg$estimate==1), level=0.5, add=TRUE, col=c("transparent", col[2]), ...))

    if (add.legend) legend(pos.legend, legend=c(expression(f[1]>f[2]), expression(f[1]<f[2])), fill=col, bty="n") 
}

plotkde.loctest.2d <- function(x, col, add=FALSE, add.legend=TRUE, pos.legend="topright", alpha=1, ...)
{ 
    if (missing(col)) col <- hcl.colors(palette="Dark2",2)
    col <- transparency.col(col, alpha=alpha)
    try(plot(x$fhat.diff.pos, col=c("transparent", col[1]), abs.cont=0.5, drawlabel=FALSE, disp="filled.contour", alpha=alpha, add=add, ...))
    try(plot(x$fhat.diff.neg, col=c("transparent", col[2]), abs.cont=0.5, drawlabel=FALSE, disp="filled.contour", alpha=alpha, add=TRUE, ...))

    if (add.legend) legend(pos.legend, legend=c(expression(f[1]>f[2]), expression(f[1]<f[2])), fill=col, bty="n")  
}

plotkde.loctest.3d <- function(x, col, color, add=FALSE, box=TRUE, axes=TRUE, alphavec=c(0.5, 0.5), add.legend=TRUE, ...)
{
    if (length(alphavec)==1) alphavec <- rep(alphavec,2)
    if (missing(col)) col <- hcl.colors(palette="Dark2",2) 
    try(plot(x$fhat.diff.pos, color=col[1], col=col[1], abs.cont=0.5, add=add, box=FALSE, axes=FALSE, alphavec=alphavec[1],  ...))
    try(plot(x$fhat.diff.neg, color=col[2], col=col[2], abs.cont=0.5, add=TRUE, box=box, axes=axes, alphavec=alphavec[2], ...)) 

    if (add.legend) 
    {
        if (!requireNamespace("plot3D", quietly=TRUE)) stop("Install the plot3D package as it is required.", call.=FALSE)
        plot3D::colkey(clim=c(1,2), col=col, at=c(1.25, 1.75), add=TRUE, side=1, addlines=TRUE, length=0.2, dist=-0.1, shift=0.1, labels=c(expression(f[1]>f[2]), expression(f[1]<f[2])))
    }
}
