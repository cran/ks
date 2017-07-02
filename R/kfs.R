###############################################################################
### Feature significance for ultivariate kernel density stimate 
###############################################################################

kfs <- function(x, H, h, deriv.order=2, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE, signif.level=0.05)
{
    r <- deriv.order
  
    if (is.vector(x))
    {
        if (missing(H)) {d <- 1; n <- length(x)}
        else
        {
            if (is.vector(H)) { d <- 1; n <- length(x)}
            else {x <- matrix(x, nrow=1); d <- ncol(x); n <- nrow(x)}
        }
    }
    else {d <- ncol(x); n <- nrow(x)}

    if (!missing(w))
        if (!(identical(all.equal(sum(w), n), TRUE)))
        {
            warning("Weights don't sum to sample size - they have been scaled accordingly\n")
      w <- w*n/sum(w)
        }
    if (missing(w)) w <- rep(1,n)
    
    if (missing(h) & d==1) h <- hpi(x=x, nstage=2, binned=TRUE, bgridsize=bgridsize, deriv.order=r)
    
    if (missing(H) & d>1) H <- Hpi(x=x, nstage=2-(d>2), binned=default.bflag(d=d, n=n), deriv.order=2, verbose=verbose)
    
    if (d==1)
    {
        fhatr <- kdde(x=x, h=h, deriv.order=r, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, deriv.vec=FALSE, verbose=verbose)

        fhat <- kde(x=x, h=h, gridsize=gridsize, gridtype=gridtype, xmin=min(fhatr$eval.points), xmax=max(fhatr$eval.points), binned=binned, bgridsize=bgridsize, positive=positive, adj.positive=adj.positive, w=w)
        
        fhat.est <- as.vector(fhat$estimate)
        fhatr.est <- as.vector(fhatr$estimate)

        RDrK <- (-1)^r*psins.1d(r=2*r, sigma=1)
        fhatr.Sigma <- n^(-1)* h^(-2*r-1)*RDrK
        fhatr.Sigma12 <- sqrt(fhatr.Sigma)
        fhatr.est <- fhatr.est/fhatr.Sigma12

        local.mode <- fhatr.est <= 0
        fhatr.wald <- fhatr.est^2

        gridsize <- length(fhat$estimate)        
    }
    else if (d>1)
    {
        ## KDE
        fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, verbose=verbose)
        
        ## KDDE for r=2
        fhatr <- kdde(x=x, H=H, deriv.order=r, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, eval.points=eval.points, binned=binned, bgridsize=bgridsize, w=w, deriv.vec=FALSE, verbose=verbose)
        
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
        
        gridsize <- dim(fhat$estimate)
    }

    ## Hochberg adjustment for sequential tests 
    pval.wald <- 1 - pchisq(fhatr.wald, d*(d+1)/2)
    pval.wald[fhat.est<=contourLevels(fhat, cont=100)] <- NA
    pval.wald.ord <- pval.wald[order(pval.wald)]
    num.test <- sum(!is.na(pval.wald.ord))
    
    if (num.test>=1) num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else num.test.seq <- rep(NA, prod(gridsize))
    reject.nonzero <- ((pval.wald.ord <= signif.level/(num.test + 1 - num.test.seq)) &(pval.wald.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)
    
    signif.wald <- array(FALSE, dim=gridsize)
    
    ## p-value == 0 => reject null hypotheses automatically
    signif.wald[which(pval.wald==0, arr.ind=TRUE)] <- TRUE
  
    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind
    for (i in reject.nonzero.ind)
        signif.wald[which(pval.wald==pval.wald.ord[i], arr.ind=TRUE)] <- TRUE 

    ## ESS = effective sample size
    ##ess <- n*fhat$estimate*dmvnorm.mixt(x=rep(0,d), mu=rep(0,d), Sigma=H, props=1)
    ##signif.ess <- ess >= 5

    signif.wald <- signif.wald & local.mode & fhat.est>contourLevels(fhat, cont=99)
    fhatr$estimate <- signif.wald+0
    ##fhatr$dens.estimate <- fhat.est
    class(fhatr) <- "kfs"
    
    return(fhatr)
}


plot.kfs <- function(x, display="filled.contour2", col="orange", colors="orange", abs.cont, alphavec=0.4, add=FALSE, ...)
{
    fhatr <- x
    fhatr$deriv.order <- NULL
    class(fhatr) <- "kde"

    if (is.vector(fhatr$H)) d <- 1 else d <- ncol(fhatr$H)
    if (d==1)
    {
        fhat <- kde(x=fhatr$x)
        plot(fhat, col="grey", add=add, ...)
        gridsize <- length(fhatr$estimate)
        sc.ind <- which(fhatr$estimate==1)
        sc.len <- length(sc.ind)
        sc.ind.diff <- diff(sc.ind)
        jump.ind <- which(sc.ind.diff!=1)
        jump.num <- length(jump.ind)
        
        if (jump.num==0) lines(fhat$eval.points[sc.ind], fhat$estimate[sc.ind],col=col, ...)
        
        if (jump.num > 0)
        {
            curr.ind <- sc.ind[1:jump.ind[1]]
            lines(fhat$eval.points[curr.ind], fhat$estimate[curr.ind], col=col, ...)
            if (jump.num > 1) 
            { 
                for (j in 2:length(jump.num))
                {
                    curr.ind <- sc.ind[(jump.ind[j-1]+1):jump.ind[j]]
                    lines(fhat$eval.points[curr.ind], fhat$estimate[curr.ind], col=col, ...)
                }
            }
            curr.ind <- sc.ind[(max(jump.ind)+1):sc.len]
            lines(fhat$eval.points[curr.ind], fhat$estimate[curr.ind], col=col, ...)
        }  
    }
    else if (d==2)
    {
        if (missing(abs.cont)) abs.cont <- 0.5
        disp1 <- match.arg(display, c("slice", "persp", "image", "filled.contour", "filled.contour2"))
        if (disp1=="filled.contour2") col <- c("transparent", col)
        if (disp1=="filled.contour")
        {
            col.fun <- function(n){return(c("transparent", rep("orange",n)))}
            plot(fhatr, abs.cont=abs.cont, drawlabels=FALSE, col.fun=col.fun, add=add, display=display, ...)
        }
        else
            plot(fhatr, abs.cont=abs.cont, drawlabels=FALSE, col=col, add=add, display=display, ...)
    }
    else if (d==3)
    {
        if (missing(abs.cont)) abs.cont <- 0.25
        plot(fhatr, abs.cont=abs.cont, colors=colors, alphavec=alphavec, add=add, ...)
    }
    invisible()
}
