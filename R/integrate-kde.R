#############################################################################
## Cumulative integral for KDE
#############################################################################

integral.kde <- function(q, fhat, density) #, exact=FALSE)
{
    gridsize <- length(fhat$eval.points)

    ## Use Simpson's rule to compute numerical integration

    simp.rule <- rep(0, gridsize-1)
    for (i in 1:(gridsize-1))
    {
        del <- fhat$eval.points[i+1] - fhat$eval.points[i]
        simp.rule[i] <- min(fhat$estimate[i], fhat$estimate[i+1])*del + 1/2*abs(fhat$estimate[i+1] - fhat$estimate[i])*del 
    }

    ## add last incomplete trapezoid 
    q.ind <- findInterval(x=q, vec=fhat$eval.points)
    q.prob <- rep(0, length(q))
    i <- 0

    for (qi in q.ind)
    {
        i <- i+1

        if (qi==0)
            q.prob[i] <- 0
        else if (qi < gridsize)
        {
            ## linearly interpolate kde 
            fhat.estqi <- (fhat$est[qi+1] - fhat$est[qi])/(fhat$eval[qi+1] - fhat$eval[qi]) * (q[i] - fhat$eval[qi]) + fhat$est[qi]
            delqi <- q[i] - fhat$eval[qi] 
              
            simp.ruleqi <- min(fhat.estqi, fhat$est[qi])*delqi + 1/2*abs(fhat.estqi - fhat$est[qi])*delqi
            q.prob[i] <- sum(simp.rule[1:qi]) + simp.ruleqi
        }
        else
        {
            if (density) q.prob[i] <- 1
            else q.prob[i] <- sum(simp.rule) 
        }
    }

    if (density) q.prob[q.prob>=1] <- 1

    return(q.prob)
}

## cumulative probability P(fhat <= q)

pkde <- function(q, fhat)
{
    return(integral.kde(q=q, fhat=fhat, density=TRUE))
}

## density value of KDE at x 
## alias for predict.kde

dkde <- function(x, fhat)
{
    return(predict(fhat, x=x))
}

## p-quantile of KDE, i.e. solve for x where P(fhat < x) = p 

qkde <- function(p, fhat)
{
    if (any(p > 1) | any(p < 0)) stop("p must be <= 1 and >= 0")
  
    cumul.prob <- pkde(q=fhat$eval.points, fhat=fhat)
    ind <- findInterval(x=p, vec=cumul.prob)

    quant <- rep(0, length(ind))
    for (j in 1:length(ind))
    {
        i <- ind[j]
        if (i==0)
            quant[j] <- fhat$eval.points[1]
        else if (i>=length(fhat$eval.points))
            quant[j] <- fhat$eval.points[length(fhat$eval.points)]
        else  
        {
            quant1 <- fhat$eval.points[i]
            quant2 <- fhat$eval.points[i+1]
              
            prob1 <- cumul.prob[i]
            prob2 <- cumul.prob[i+1]
            alpha <- (p[j] - prob2)/(prob1 - prob2)
            quant[j] <- quant1*alpha + quant2*(1-alpha)
        }
    }
  
    return(quant)
}

## Silverman (1983)'s random sample from KDE

rkde <- function(n, fhat, positive=FALSE)
{
    if (positive) x <- log(fhat$x)
    else x <- fhat$x

    if (is.vector(fhat$H)) { d <- 1; nsamp <- length(x) }
    else { d <- ncol(fhat$H); nsamp <- nrow(x) } 

    x.ind <- sample(1:nsamp, size=n, replace=TRUE, prob=fhat$w)

    if (d==1)
    {
        h <- fhat$h
        rkde.val <- x[x.ind] + h*rnorm(n=n, mean=0, sd=1)
        if (positive) rkde.val <- exp(rkde.val)
    }
    else if (d>1)
    {
        H <- fhat$H
        rkde.val <- x[x.ind,] + rmvnorm(n=n, mean=rep(0,d), sigma=diag(d)) %*% matrix.sqrt(H)
    }

    return(rkde.val)
}

## plot cumulative probability as shaded region on a KDE 

plotkde.cumul <- function(fhat, q, add=FALSE, col="blue", ...)
{
    qind <- fhat$eval.points<=q
    n <- sum(qind)
    if (!add) plot(fhat)
    polygon(c(fhat$eval.points[qind],fhat$eval.points[n]), c(fhat$estimate[qind],0), col=col, ...)
}
