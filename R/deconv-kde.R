######################################################################
## Deconvolution KDE
######################################################################
dckde <- function(...) { return (kdcde(...)) }

kdcde <- function(x, H, h, Sigma, sigma, reg, bgridsize, gridsize, binned, verbose=FALSE, ...)
{
    ## default values 
    ksd <- ks.defaults(x=x, binned=binned, bgridsize=bgridsize, gridsize=gridsize)
    d <- ksd$d; n <- ksd$n 
    binned <- ksd$binned
    gridsize <- ksd$gridsize
    bgridsize <- ksd$bgridsize

    x <- as.matrix(x)
    if (d==1 & missing(h)) h <- hpi(x=x, nstage=2, binned=binned, deriv.order=0)
    if (d>1 & missing(H)) H <- Hpi(x=x, nstage=2, binned=binned, deriv.order=0)

    if (d==1) stop("dckde not yet implemented for d=1")
    if (missing(reg)) reg <- reg.ucv(x=x, H=H, h=h, Sigma=Sigma, sigma=sigma, k=5, d=d, binned=binned, verbose=verbose)

    ## Deconvolution KDE is weighted KDE with non-uniform weights
    w <- dckde.weights(x=x, H=H, Sigma=Sigma, reg=reg)
    fhat <- kde(x=x, H=H, w=w, binned=binned, bgridsize=bgridsize, gridsize=gridsize, verbose=verbose, ...)
    fhat$reg <- reg

    return(fhat)
}

## weights for deconvolution KDE
## code adapted from DeconWK 0.6-5, author B. Turlach
## R-forge website: https://r-forge.r-project.org/R/?group_id=630
dckde.weights <- function(x, Sigma, H, reg)
{
    n <- nrow(x)
    d <- ncol(x)
    Qmat <- matrix(0, ncol=n, nrow=n)
    bvec <- rep(0, n)
    
    for (j in 1:n)
    {  
        Qmat[j,] <- dmvnorm.mixt(x, mus=x[j,], Sigmas=2*H + 2*Sigma, props=1)
        bvec[j] <- sum(dmvnorm.mixt(x, mus=x[j,], Sigmas=2*H + Sigma, props=1))
    }
    
    if(!missing(reg)) diag(Qmat) <- diag(Qmat) + reg/n
    bvec <- bvec/n
    
    val <- kernlab::ipop(c=-bvec, H=Qmat, A=rep(1,n), b=1, r=0, l=rep(0,n), u=rep(1, n)) 
    w <- kernlab::primal(val)
    w <- w/sum(w)*n
    return(w)  
}

## unbiased k-fold cross validation choice of regularisation penalty (gamma)
reg.ucv <- function(x, H, h, Sigma, sigma, k=5, d, binned=FALSE, verbose=FALSE)
{
    gamma.ucv.temp <- function(gamma) { return(-reg.ucv.val(x=x, H=H, Sigma=Sigma, k=k, reg=gamma^2, binned=binned)) }
    gamma.val <- nlm(f=gamma.ucv.temp, p=0.1, print.level=2*as.numeric(verbose))$estimate^2
    
    return(gamma.val)
}

## k-fold UCV value for regularisation penalty
reg.ucv.val <- function(x, Sigma, H, reg, k=5, binned=FALSE)
{
    n <- nrow(x)
    n.seq <- block.indices(n, n, npergroup=round(n/k))
    cv.val <- 0
    for (j in 1:(length(n.seq)-1))
    {
        iind <- n.seq[j]:(n.seq[j+1]-1)
        w <- dckde.weights(x=x[-iind,], H=H, Sigma=Sigma, reg=reg)  
        fhat <- kde(x=x[-iind,], H=H+Sigma, w=w, binned=binned)
        cv.val <- cv.val + sum(predict(fhat, x=x[iind,]))
    }
    return(cv.val)
}
