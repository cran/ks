#####################################################################
## Kernel estimators of the multivariate cdf (cumulative distribution function)
#####################################################################

kcde <- function(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE, tail.flag="lower.tail")
{
  r <- 0
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
  {
    if (!(identical(all.equal(sum(w), n), TRUE)))
    {
      warning("Weights don't sum to sample size - they have been scaled accordingly\n")
      w <- w*n/sum(w)
    }
  }
  else w <- rep(1,n)

  tail.flag1 <- match.arg(tail.flag, c("lower.tail", "upper.tail")) 

  ## KCDE is computed as cumulative Riemann sum of KDE on a grid
  if (d==1)
  {
    if (missing(h)) h <- hpi.kcde(x=x, binned=default.bflag(d=d, n=n))
    Fhat <- kde(x=x, h=h, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, binned=binned, bgridsize=bgridsize, positive=positive, adj.positive=adj.positive, w=w)
    diffe <- abs(diff(Fhat$eval.points))
    if (tail.flag1=="lower.tail") Fhat$estimate <- c(0, diffe) * cumsum(Fhat$estimate)
    else Fhat$estimate <- c(diffe[1], diffe) * (sum(Fhat$estimate) - cumsum(Fhat$estimate))
  }
  else if (d==2)
  {
    if (missing(H)) Hpi.kcde(x=x, binned=default.bflag(d=d, n=n), bgridsize=bgridsize, verbose=FALSE)
    Fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, binned=binned, bgridsize=bgridsize, w=w)
    diffe1 <- abs(diff(Fhat$eval.points[[1]]))
    diffe2 <- abs(diff(Fhat$eval.points[[2]]))

    if (tail.flag1=="lower.tail")
    {
      Fhat$estimate <- apply(Fhat$estimate, 1, cumsum)*c(0,diffe1)
      Fhat$estimate <- apply(t(Fhat$estimate), 2, cumsum)*c(0,diffe2)
    }
    else
    {
      Fhatsum <- matrix(apply(Fhat$estimate, 1, sum), ncol=ncol(Fhat$estimate), nrow=nrow(Fhat$estimate), byrow=TRUE)
      Fhat$estimate <- (Fhatsum-apply(Fhat$estimate, 1, cumsum))*c(diffe1[1], diffe1)
      Fhatsum <- matrix(apply(Fhat$estimate, 1, sum), ncol=ncol(Fhat$estimate), nrow=nrow(Fhat$estimate), byrow=TRUE)
      Fhat$estimate <- (Fhatsum-apply(t(Fhat$estimate), 2, cumsum))*c(diffe2[1], diffe2)
    }
  }
  else if (d==3)
  {
     if (missing(H)) Hpi.kcde(x=x, binned=default.bflag(d=d, n=n), bgridsize=bgridsize, verbose=FALSE)
     Fhat <- kde(x=x, H=H, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, binned=binned, bgridsize=bgridsize, w=w)
    Fhat.temp <- Fhat$estimate
    diffe1 <- abs(diff(Fhat$eval.points[[1]]))
    diffe2 <- abs(diff(Fhat$eval.points[[2]]))
    diffe3 <- abs(diff(Fhat$eval.points[[3]]))
    if (tail.flag1=="lower.tail")
    {
      for (i in 1:dim(Fhat$estimate)[3])
      {
        Fhat.temp[,,i] <- apply(Fhat.temp[,,i], 1, cumsum)*c(0,diffe1)
        Fhat.temp[,,i] <- apply(t(Fhat.temp[,,i]), 2, cumsum)*c(0,diffe2)
      }
      for (i in 1:dim(Fhat$estimate)[1])
        for (j in 1:dim(Fhat$estimate)[2])
          Fhat.temp[i,j,] <- cumsum(Fhat.temp[i,j,])*c(0,diffe3)
      Fhat$estimate <- Fhat.temp
    }
    else
    {
      for (i in 1:dim(Fhat$estimate)[3])
      {
         Fhatsum <- matrix(apply(Fhat.temp[,,i], 1, sum), ncol=ncol(Fhat.temp), nrow=nrow(Fhat.temp), byrow=TRUE)
         Fhat.temp[,,i] <- (Fhatsum-apply(Fhat.temp[,,i], 1, cumsum))*c(diffe1[1], diffe1)
         Fhatsum <- matrix(apply(Fhat.temp[,,i], 1, sum), ncol=ncol(Fhat.temp), nrow=nrow(Fhat.temp), byrow=TRUE)
        Fhat.temp[,,i] <- (Fhatsum-apply(t(Fhat.temp[,,i]), 2, cumsum))*c(diffe2[1],diffe2)
      }
      
      for (i in 1:dim(Fhat$estimate)[1])
        for (j in 1:dim(Fhat$estimate)[2])
          {
            Fhatsum <- sum(Fhat.temp[i,j,])
            Fhat.temp[i,j,] <- (Fhatsum-cumsum(Fhat.temp[i,j,]))*c(diffe3[1],diffe3)
          }
      Fhat$estimate <- Fhat.temp
    }
  }
  ## normalise max CDF estimate equal to 1
  Fhat$estimate <- Fhat$estimate/max(Fhat$estimate)

  if (!missing(eval.points))
  {
    if (d<=3)
    {
      Fhat$estimate <- predict(Fhat, x=eval.points)
      Fhat$eval.points <- eval.points
    }
    else
    {
      Fhat <- kcde.points(x=x, H=H, eval.points=eval.points, w=w, verbose=verbose, tail.flag=tail.flag1)
    }
  }
  Fhat$tail <- tail.flag1
  class(Fhat) <- "kcde"
  return(Fhat) 
}


## KCDE is computed at specified estimation points

kcde.points <- function(x, H, eval.points, w, verbose=FALSE, tail.flag="lower.tail") 
{
  n <- nrow(x)
  if (missing(w)) w <- rep(1,n)
  if (verbose) pb <- txtProgressBar() 
  Fhat <- rep(0, nrow(eval.points))
  pmvnorm.temp <- function(x, ...) { return(pmvnorm(mean=x, ...)) }
   
  for (i in 1:nrow(eval.points))
  {  
    if (verbose) setTxtProgressBar(pb, i/(nrow(eval.points)-1))
    if (tail.flag=="lower.tail")
      Fhat[i] <- sum(apply(x, 1, pmvnorm.temp, upper=eval.points[i,], sigma=H))
    else
      Fhat[i] <- sum(apply(x, 1, pmvnorm.temp, lower=eval.points[i,], sigma=H))
  }
  Fhat <- Fhat/n
  if (verbose) close(pb)
  
  return(list(x=x, eval.points=eval.points, estimate=Fhat, H=H, gridded=FALSE, binned=FALSE, names=NULL, w=w))
}



#####################################################################
## Plotting functions for 1-d to 3-d KCDE
#####################################################################

plot.kcde <- function(x, ...)
{ 
  Fhat <- x
  if (is.vector(Fhat$x)) plotkcde.1d(Fhat, ...)
  else
  {
    d <- ncol(Fhat$x)

    if (d==2) 
    {
      plotret <- plotkcde.2d(Fhat, ...)
      invisible(plotret)
    }
    else if (d==3)
    {
      plotkcde.3d(Fhat, ...)
      invisible()
    }
    else stop ("plot function only available for 1, 2 or 3-d data")
  }
}


plotkcde.1d <- function(Fhat, xlab, ylab="Distribution function", add=FALSE, drawpoints=FALSE, ptcol="blue", jitter=FALSE, ...) 
{
  if (missing(xlab)) xlab <- Fhat$names
  if (Fhat$tail=="upper.tail") zlab <- "Survival function"
  if (add) lines(Fhat$eval.points, Fhat$estimate, xlab=xlab, ylab=ylab, ...)
  else plot(Fhat$eval.points, Fhat$estimate, type="l", xlab=xlab, ylab=ylab, ...) 

  if (drawpoints)
    if (jitter) rug(jitter(Fhat$x), col=ptcol)
    else rug(Fhat$x, col=ptcol)
}


plotkcde.2d <- function(Fhat, display="persp", cont=seq(10,90, by=10), abs.cont,
    xlab, ylab, zlab="Distribution function", cex=1, pch=1, labcex,  
    add=FALSE, drawpoints=FALSE, drawlabels=TRUE, theta=-30, phi=40, d=4,
    ptcol="blue", col, lwd=1, border=NA, thin=1, ...) 
{
  disp1 <- match.arg(display, c("slice", "persp", "image", "filled.contour", "filled.contour2"))
  
  if (!is.list(Fhat$eval.points)) stop("need a grid of density estimates")

  if (missing(xlab)) xlab <- Fhat$names[1]
  if (missing(ylab)) ylab <- Fhat$names[2]
  if (missing(labcex)) labcex <-1
  if (Fhat$tail=="upper.tail") zlab <- "Survival function"
  
  ## perspective/wireframe plot
  if (disp1=="persp")
  {
    hts <- seq(0, 1.1*max(Fhat$estimate), length=100)
    if (missing(col)) col <- grey(seq(0,0.9, length=length(hts)+1)) ## rev(heat.colors(length(hts)+1)) #
    if (length(col)<100) col <- rep(col, length=100)

    ## thinning indices
    plot.ind <- list(seq(1, length(Fhat$eval.points[[1]]), by=thin), seq(1, length(Fhat$eval.points[[2]]), by=thin))

    z <- Fhat$estimate[plot.ind[[1]], plot.ind[[2]]]
    nrz <- nrow(z)
    ncz <- ncol(z)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, length(hts)+1)
    
    plotret <- persp(Fhat$eval.points[[1]][plot.ind[[1]]], Fhat$eval.points[[2]][plot.ind[[2]]], z, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, col=col[facetcol], border=border, ...)
  }
  else if (disp1=="slice") 
  {
    if (!add)
      plot(Fhat$x[,1], Fhat$x[,2], type="n", xlab=xlab, ylab=ylab, ...)
 
    if (missing(abs.cont)) hts <- cont/100
    else hts <- abs.cont
    
    if (missing(col)) col <- 1 
    if (length(col)<length(hts)) col <- rep(col, times=length(hts))
    
    ## draw contours         
    for (i in 1:length(hts)) 
    {
      if (missing(abs.cont)) scale <- cont[i]/hts[i]
      else scale <- 1
      if (hts[i]>0)
        contour(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate*scale, level=hts[i]*scale, add=TRUE, drawlabels=drawlabels, labcex=labcex, col=col[i], lwd=lwd, ...)
    }
    ## add points 
    if (drawpoints) points(Fhat$x[,1], Fhat$x[,2], col=ptcol, cex=cex, pch=pch)
  }
  ## image plot
  else if (disp1=="image")
  {
    image(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate, xlab=xlab, ylab=ylab, add=add, ...)
    box()
  }
  else if (disp1=="filled.contour" | disp1=="filled.contour2") 
  {
    hts <- cont/100
    if (missing(col)) col <- c("transparent", rev(heat.colors(length(hts))))
    
    clev <- c(-0.01*max(abs(Fhat$estimate)), hts, max(c(Fhat$estimate, hts)) + 0.01*max(abs(Fhat$estimate)))
    
    if (disp1=="filled.contour2")
    {
      image(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate, xlab=xlab, ylab=ylab, add=add, col=col[1:(length(hts)+1)], breaks=clev, ...)

      ## draw contours         
      for (i in 1:length(hts)) 
        contour(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate, level=hts[i], add=TRUE, drawlabels=FALSE, col=col[i+1], lwd=7)
      if (!missing(lwd))
      {
        for (i in 1:length(hts)) 
        {
          if (missing(abs.cont)) scale <- cont[i]/hts[i]
          else scale <- 1
          
          if (lwd >=1) contour(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate*scale, level=hts[i]*scale, add=TRUE, drawlabels=drawlabels, col=1, labcex=labcex, lwd=lwd, ...)
        }
      }
      ## add points 
      if (drawpoints) points(Fhat$x[,1], Fhat$x[,2], col=ptcol, cex=cex, pch=pch)
    }
    else
    {
      if (tail(hts, n=1) < max(Fhat$estimate)) hts <- c(hts, max(Fhat$estimate))
      filled.contour(Fhat$eval.points[[1]], Fhat$eval.points[[2]], Fhat$estimate, xlab=xlab, ylab=ylab, levels=hts, ...)
    }
  }
  if (disp1=="persp")  invisible(plotret)
  else invisible()
}
  

plotkcde.3d <- function(Fhat, cont=c(25,50,75), colors, alphavec, size=3, ptcol="blue", add=FALSE, xlab, ylab, zlab, drawpoints=FALSE, alpha=1, box=TRUE, axes=TRUE, ...)

{
  hts <- sort(cont/100)
  nc <- length(hts)
  
  if (missing(colors)) colors <- rev(heat.colors(nc))
  if (missing(xlab)) xlab <- Fhat$names[1]
  if (missing(ylab)) ylab <- Fhat$names[2]
  if (missing(zlab)) zlab <- Fhat$names[3]
  if (missing(alphavec)) alphavec <- seq(0.5,0.1,length=nc)

  if (drawpoints)
    plot3d(Fhat$x[,1],Fhat$x[,2],Fhat$x[,3], size=size, col=ptcol, alpha=alpha, xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
  else
    plot3d(Fhat$x[,1],Fhat$x[,2],Fhat$x[,3], type="n", xlab=xlab, ylab=ylab, zlab=zlab, add=add, box=FALSE, axes=FALSE, ...)
  bg3d(col="white")
  
  for (i in 1:nc)
    if (hts[nc-i+1] < max(Fhat$estimate))
      contour3d(Fhat$estimate, level=hts[nc-i+1], x=Fhat$eval.points[[1]], y=Fhat$eval.points[[2]], z=Fhat$eval.points[[3]], add=TRUE, color=colors[i], alpha=alphavec[i], box=FALSE, axes=FALSE, ...)

  if (box) box3d()
  if (axes) axes3d()
}



#####################################################################
## Bandwidth selectors for KCDE
#####################################################################

### Normal scale bandwidth selectors

hns.kcde <- function(x)
{
  d <- 1
  n <- length(x)
  #m1 <- 0.2820948
  sigma <- sd(x)
  hns <- 4^(1/3)*sigma*n^(-1/3)

  return(hns)
}

Hns.kcde <- function(x)
{
  if (is.vector(x)) {return(hns.kcde(x)^2)}
  d <- ncol(x)
  n <- nrow(x)
  m1 <- (4*pi)^(-1/2)

  Jd <- matrix(1, ncol=d, nrow=d)
  Sigma <- var(x)
  Hns <- (4*det(Sigma)^(1/2)*tr(matrix.sqrt(Sigma))/tr(Sigma))^(2/3)*Sigma*n^(-2/3)

  return(Hns)
}



## Plug-in bandwidth selector

hpi.kcde <- function(x, nstage=2, binned=TRUE)
{
  n <- length(x)
  d <- 1
  
  K2 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=2)  
  K4 <- dnorm.deriv(x=0, mu=0, sigma=1, deriv.order=4) 
  m2 <- 1  
  m1 <- (4*pi)^(-1/2)
  
  ## formula for bias annihliating bandwidths from Wand & Jones (1995, p.70)
  if (nstage==2)
  {
    psi6.hat <- psins.1d(r=6, sigma=sd(x))
    gamse4 <- (2*K4/(-m2*psi6.hat*n))^(1/(4+3)) 
    psi4.hat <- kfe.1d(x=x, g=gamse4, deriv.order=4, inc=1, binned=binned)
    gamse2 <- (2*K2/(-m2*psi4.hat*n))^(1/(2+3))
    psi2.hat <- kfe.1d(x=x, g=gamse2, deriv.order=2, inc=1, binned=binned)
  }
  else 
  {
    psi4.hat <- psins.1d(r=4, sigma=sd(x))
    gamse2 <- (2*K2/(-m2*psi4.hat*n))^(1/(2+3))
    psi2.hat <- kfe.1d(x=x, g=gamse2, deriv.order=2, inc=1, binned=binned)
  }

  ## formula form Polansky & Baker (2000)
  h <- (2*m1/(-m2^2*psi2.hat*n))^(1/3) 

  return(h)
}


Hpi.kcde <- function(x, nstage=2, pilot, Hstart, binned=FALSE, bgridsize, amise=FALSE, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  m1 <- (4*pi)^(-1/2)
  Jd <- matrix(1, ncol=d, nrow=d)
  
  if(!is.matrix(x)) x <- as.matrix(x)
  if (missing(pilot)) pilot <- "dunconstr"
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))
 
  if (pilot1=="dscalar") stop("use dunconstr pilot for Hpi.kcde since pre-scaling approaches are not valid")
  
  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=4, Sigma=var(x), deriv.vec=TRUE)
    
    amse2.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      Hinv <- chol2inv(chol(H))
      Hinv12 <- matrix.sqrt(Hinv)
      amse2.val <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
      return(sum(amse2.val^2)) 
    }
      
    Hstart2 <- matrix.sqrt(Gns(r=2, n=n, Sigma=var(x)))
    optim.fun1 <- match.arg(optim.fun, c("nlm", "optim"))
 
    if (optim.fun1=="nlm")
    {
      result <- nlm(p=vech(Hstart2), f=amse2.temp, print.level=2*as.numeric(verbose))    
      H2 <- invvech(result$estimate) %*% invvech(result$estimate)
    }
    else
    {
      result <- optim(vech(Hstart2), amse2.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
      H2 <- invvech(result$par) %*% invvech(result$par)
    }
 
    psi2.hat <- kfe(x=x, G=H2, deriv.order=2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }
  else
    psi2.hat <- psins(r=2, Sigma=var(x), deriv.vec=TRUE)    
  
  if (missing(Hstart)) Hstart <- Hns.kcde(x=x)
  
  ## stage 2
  amise.temp <- function(vechH)
  { 
    H <- invvech(vechH) %*% invvech(vechH)
    H12 <- matrix.sqrt(H)
    amise.val <- -2*n^(-1)*m1*tr(H12) - 1/4*t(vec(H %*% H)) %*% psi2.hat
    ##amise.val <- -2*n^(-1)*m1*sum(H12) - 1/4*t(vec(H %*% H)) %*% psi2.hat
    return(drop(amise.val)) 
  }
  
  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- match.arg(optim.fun, c("optim", "nlm"))
  
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=vech(Hstart), f=amise.temp, print.level=2*as.numeric(verbose)) 
    H <- invvech(result$estimate) %*% invvech(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(vech(Hstart), amise.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- invvech(result$par) %*% invvech(result$par)
    amise.star <- result$value
  }

  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}

Hpi.diag.kcde <- function(x, nstage=2, pilot, Hstart, binned=FALSE, bgridsize, amise=FALSE, verbose=FALSE, optim.fun="nlm")
{
  n <- nrow(x)
  d <- ncol(x)
  m1 <- (4*pi)^(-1/2)
  Jd <- matrix(1, ncol=d, nrow=d)
  
  if(!is.matrix(x)) x <- as.matrix(x)
  if (missing(pilot)) pilot <- "dscalar"
  pilot1 <- match.arg(pilot, c("dunconstr", "dscalar"))
  if (pilot1=="dunconstr") stop("use dscalar pilot for Hpi.diag.kcde since pre-sphering approaches are not valid")

  D2K0 <- t(dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=2))
  if (nstage==2)
  {  
    ## stage 1
    psi4.ns <- psins(r=4, Sigma=var(x), deriv.vec=TRUE)
    
    amse2.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      Hinv <- chol2inv(chol(H))
      Hinv12 <- matrix.sqrt(Hinv)
      amse2.val <- 1/(det(H)^(1/2)*n)*((Hinv12 %x% Hinv12) %*% D2K0) + 1/2* t(vec(H) %x% diag(d^2)) %*% psi4.ns
      return(sum((amse2.val)^2)) 
    }
      
    Hstart2 <- matrix.sqrt(Gns(r=2, n=n, Sigma=var(x)))
    optim.fun1 <- match.arg(optim.fun, c("optim", "nlm")) 
 
    if (optim.fun1=="nlm")
    {
      result <- nlm(p=diag(Hstart2), f=amse2.temp, print.level=2*as.numeric(verbose))    
      H2 <- diag(result$estimate) %*% diag(result$estimate)
    }
    else
    {
      result <- optim(diag(Hstart2), amse2.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
      H2 <- diag(result$par) %*% diag(result$par)
    }
 
    psi2.hat <- kfe(x=x, G=H2, deriv.order=2, add.index=FALSE, binned=binned, bgridsize=bgridsize, verbose=verbose)
  }
  else
    psi2.hat <- psins(r=2, Sigma=var(x), deriv.vec=TRUE)    
  
  if (missing(Hstart)) Hstart <- Hns.kcde(x=x)
  
  ## stage 2
  amise.temp <- function(diagH)
  { 
    H <- diag(diagH) %*% diag(diagH)
    H12 <- matrix.sqrt(H)
    amise.val <- -2*n^(-1)*m1*tr(H12) - 1/4*t(vec(H %*% H)) %*% psi2.hat
    return(drop(amise.val)) 
  }
  
  Hstart <- matrix.sqrt(Hstart)
  optim.fun1 <- match.arg(optim.fun, c("optim", "nlm")) 
  
  if (optim.fun1=="nlm")
  {
    result <- nlm(p=diag(Hstart), f=amise.temp, print.level=2*as.numeric(verbose)) 
    H <- diag(result$estimate) %*% diag(result$estimate)
    amise.star <- result$minimum
  }
  else
  {
    result <- optim(diag(Hstart), amise.temp, method="BFGS", control=list(trace=as.numeric(verbose)))
    H <- diag(result$par) %*% diag(result$par)
    amise.star <- result$value
  }

  if (!amise) return(H)
  else return(list(H=H, PI=amise.star))
}

#####################################################################
## Multivariate kernel ROC estimators
#####################################################################

kroc <- function(x1, x2, H1, h1, hy, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE, nref=1e4)
{
  if (is.vector(x1)) {d <- 1; n1 <- length(x1)} else {d <- ncol(x1); n1 <- nrow(x1)}
  if (!missing(eval.points)) stop("eval.points in kroc not yet implemented")
  
  if (d==1)
  {
    if (missing(h1)) h1 <- hpi.kcde(x=x1, binned=default.bflag(d=d, n=n1))
    
    y <- kcde(x=x1, h=h1, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, binned=binned, bgridsize=bgridsize, positive=positive, adj.positive=adj.positive, eval.points=x2, w=w, tail.flag="upper.tail")$estimate
  }
  else
  {
    if (missing(H1)) H1 <- Hpi.kcde(x=x1, binned=default.bflag(d=d, n=n1))
    y <- kcde(x=x1, H=H1, gridsize=gridsize, gridtype=gridtype, xmin=xmin, xmax=xmax, supp=supp, binned=binned, bgridsize=bgridsize, eval.points=x2, w=w, tail.flag="upper.tail")$estimate
  }

  ## transform from [0,1] to reals
  y <- qnorm(y[y>0])
  hy <- hpi.kcde(y, binned=default.bflag(d=d, n=n1))
  Fhaty <- kcde(x=y, h=hy, binned=TRUE) 
  Fhaty$eval.points <- pnorm(Fhaty$eval.points)
  Fhaty$x <- list(x1, x2)

 
  if (d==1) {Fhaty$h1 <- h1; Fhaty$H1 <- h1^2; Fhaty$hy <- hy}
  else {Fhaty$H1 <- H1; Fhaty$hy <- hy}
  
  transform <- TRUE
  if (transform & d>1) Fhaty <- kroc.transform(Fhaty, nref=nref)
  else
  {
    ## Use spline to smooth out transformed ROC curve
    Fhaty.smoothed <- smooth.spline(Fhaty$eval.points, Fhaty$estimate)
    Fhaty.smoothed <- predict(Fhaty.smoothed, x=seq(0,1,length=length(Fhaty$eval.points))) 
    Fhaty$eval.points <- Fhaty.smoothed$x
    Fhaty$estimate <- Fhaty.smoothed$y

    ## add (0,0) and (1,1) as endpoints
    if (head(Fhaty$eval.points, n=1)!=0) Fhaty$eval.points[1] <- 0
    if (head(Fhaty$estimate, n=1)!=0) Fhaty$estimate[1] <- 0
    if (tail(Fhaty$eval.points, n=1)!=1) Fhaty$eval.points[length(Fhaty$eval.points)] <- 1
    if (tail(Fhaty$estimate, n=1)!=1) Fhaty$estimate[length(Fhaty$estimate)] <- 1
    Fhaty$estimate[Fhaty$estimate>1] <- 1
    Fhaty$estimate[Fhaty$estimate<0] <- 0
    Fhaty$indices <- indices.kroc(Fhaty)
  }
  
  Fhaty <- Fhaty[-c(4,5)] 
  class(Fhaty) <- "kroc"
  return(Fhaty)
}


## ROC curve to compare normal r.v. with variance Sigma against itself
roc.ns <- function(Sigma, nref=1e4)
{
  if (is.vector(Sigma)) d <- 1 else d <- ncol(Sigma)
  ## locally set random seed not to interfere with global random number generators
  if (!exists(".Random.seed")) rnorm(1)
  old.seed <- .Random.seed
  on.exit( { .Random.seed <<- old.seed } )
  set.seed(8192)
  xref <- rmvnorm.mixt(n=nref, mus=rep(0,d), Sigmas=Sigma)

  ##approx.ref1 <- match.arg(approx.ref, c("kcde")) 
  ##if (approx.ref1=="pmvnorm")
  ##{
  ##  y <- rep(0, nref)
  ##  for (i in 1:length(y)) y[i] <- pmvnorm(lower=xref[i,], sigma=Sigma)
  ##}
  ##else if (approx.ref1=="kcde")
  ##{
    ##y <- kcde(x=xref, binned=TRUE, H=diag(diag(Hns(xref))), tail.flag="upper.tail", eval.points=xref)$estimate
  y <- kcde(x=xref, binned=FALSE, H=Hns.kcde(xref), tail.flag="upper.tail", eval.points=xref)$estimate
  ##}
 
  ## transform from [0,1] to reals
  y <- qnorm(y[y>0])
  hy <- hpi.kcde(y, binned=default.bflag(d=d, n=length(y)))
  Fhaty <- kcde(x=y, h=hy, binned=TRUE)
  Fhaty$eval.points <- pnorm(Fhaty$eval.points)
  
  ## Use spline to smooth out transformed ROC curve
  Fhaty.smoothed <- smooth.spline(Fhaty$eval.points, Fhaty$estimate)
  Fhaty.smoothed <- predict(Fhaty.smoothed, x=seq(0,1,length=length(Fhaty$eval.points)))
  Fhaty$eval.points <- Fhaty.smoothed$x
  Fhaty$estimate <- Fhaty.smoothed$y

  ## add (0,0) and (1,1) as endpoints
  if (head(Fhaty$eval.points, n=1)!=0) Fhaty$eval.points[1] <- 0
  if (head(Fhaty$estimate, n=1)!=0) Fhaty$estimate[1] <- 0
  if (tail(Fhaty$eval.points, n=1)!=1) Fhaty$eval.points[length(Fhaty$eval.points)] <- 1
  if (tail(Fhaty$estimate, n=1)!=1) Fhaty$estimate[length(Fhaty$estimate)] <- 1

  Fhaty$estimate[Fhaty$estimate>1] <- 1
  Fhaty$estimate[Fhaty$estimate<0] <- 0
  
  Fhaty$hy <- hy
  Fhaty$indices <- indices.kroc(Fhaty)
  class(Fhaty) <- "kroc"
   
  return(Fhaty)
}


kroc.transform <- function(Rhat, Rhat.ref, nref)
{
  ## transform given ROC curve to reference ROC curve
  Rhat.trans <- Rhat
  
  if (missing(Rhat.ref))
  {
    xref <- rbind(Rhat$x[[1]], Rhat$x[[2]]) 
    Rhat.ref <- roc.ns(Sigma=var(xref), nref=nref)
  }

  Rhat.ref$estimate <- predict(Rhat.ref, x=Rhat$eval.points)
  Rhat.trans$eval.points <- Rhat.ref$estimate

  ## Use spline to smooth out transformed ROC curve
  Rhat.trans.smoothed <- smooth.spline(Rhat.trans$eval.points, Rhat.trans$estimate)
  Rhat.trans.smoothed <- predict(Rhat.trans.smoothed, x=seq(0,1,length=length(Rhat$eval.points)))
  Rhat.trans$eval.points <- Rhat.trans.smoothed$x
  Rhat.trans$estimate <- Rhat.trans.smoothed$y

  ## add (0,0) and (1,1) as endpoints
  if (head(Rhat.trans$eval.points, n=1)!=0) Rhat.trans$eval.points[1] <- 0
  if (head(Rhat.trans$estimate, n=1)!=0) Rhat.trans$estimate[1] <- 0
  if (tail(Rhat.trans$eval.points, n=1)!=1) Rhat.trans$eval.points[length(Rhat.trans$eval.points)] <- 1
  if (tail(Rhat.trans$estimate, n=1)!=1) Rhat.trans$estimate[length(Rhat.trans$estimate)] <- 1 
  Rhat.trans$estimate[Rhat.trans$estimate>1] <- 1
  Rhat.trans$estimate[Rhat.trans$estimate<0] <- 0
  
  Rhat.trans$indices <- indices.kroc(Rhat.trans)
  
  return(Rhat.trans)
}


### summary measure of ROC curves

indices.kroc <- function(Rhat)
{
  auc <- sum(abs((head(Rhat$estimate, n=-1) - tail(Rhat$estimate, n=-1)))*abs(diff(Rhat$eval.points))/2 + head(Rhat$estimate, n=-1)*abs(diff(Rhat$eval.points)))

  youden.val <- Rhat$estimate - Rhat$eval.points
  if (max(youden.val)>0.001)
  {  
    youden.ind <- which.max(youden.val)
    youden <- youden.val[youden.ind]
    LR <- list(minus=(1-Rhat$estimate[youden.ind])/(1-Rhat$eval.points[youden.ind]), plus=Rhat$estimate[youden.ind]/Rhat$eval.points[youden.ind])
  }
  else
    LR <- list(minus=1, plus=1)
  return(list(auc=auc, youden=max(youden.val), LR=LR))
}  



## plot method
plot.kroc <- function(x, add=FALSE, add.roc.ref=FALSE, ylab="True positive rate (sensitivity)", xlab=
  expression("False positive rate"~~group("(", list(bar(specificity)), ")")), ...)
{
  Rhat <- x
  if (add) lines(Rhat$eval.points, Rhat$estimate, ...)
  else plot(Rhat$eval.points, Rhat$estimate, type="l", ylab=ylab, xlab=xlab, ...)
  
  if (is.vector(Rhat$x[[1]])) d <- 1 else d <- ncol(Rhat$x[[1]])
  if (add.roc.ref)
  {
    if (d==1)
    {
      z <- seq(0,1, length=401)
      kind <- 0:(d-1)
      roc.indep <- 0
      for (k in kind) roc.indep <- roc.indep + z*(-log(z))^k/factorial(k)
      lines(z, roc.indep, lty=2, col="grey")
    }
    else
    {  
      xref <- rbind(Rhat$x[[1]], Rhat$x[[2]])
      Rhat.ref <- roc.ns(Sigma=var(xref))
      lines(Rhat.ref$eval.points, Rhat.ref$estimate, lty=2, col="grey")
    }
  }
}

## summary method
summary.kroc <- function(object, ...)
{
  cat("Summary measures for ROC curve\nAUC =", signif(object$indices$auc, ...), "\n")
  cat("Youden index =", signif(object$indices$youden, ...), "\n")
  cat(paste("(LR-, LR+) = (",  signif(object$indices$LR$minus, ...), ", ", signif(object$indices$LR$plus, ...),")\n\n",sep=""))
}

