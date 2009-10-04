
###############################################################################
## Univariate mixture normal densities
###############################################################################


rnorm.mixt <- function(n=100, mus=0, sigmas=1, props=1, mixt.label=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (mixt.label)
      rand <- cbind(rnorm(n=n, mean=mus, sd=sigmas), rep(1, n))
    else
      rand <- rnorm(n=n, mean=mus, sd=sigmas)
  }
  ### multiple component mixture
  else
  {
    k <- length(props)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## alternative method for component membership
    ##runif.memb <- runif(n=n)
    ##memb <- findInterval(runif.memb, c(0,cumsum(props)), rightmost.closed=TRUE)
    ##n.prop <- table(memb)
    
    ## compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k) ##for (i in as.numeric(rownames(n.prop)))
    {
      ## compute random sample from normal mixture component
      if (n.prop[i] > 0)
        if (mixt.label)
          rand <- rbind(rand, cbind(rnorm(n=n.prop[i], mean=mus[i], sd=sigmas[i]), rep(i, n.prop[i])))
        else
          rand <- c(rand, rnorm(n=n.prop[i], mean=mus[i], sd=sigmas[i]))
    }
  }
  if (mixt.label)
    return(rand[sample(n),])
  else
    return(rand[sample(n)])
}


dnorm.mixt <- function(x, mus=0, sigmas=1, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm(x, mean=mus[1], sd=sigmas[1])

  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0

    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dnorm(x, mean=mus[i], sd=sigmas[i])
  }
  
  return(dens)
}   


###############################################################################
## Partial derivatives of the univariate normal (mean 0) 
## 
## Parameters
## x - points to evaluate at
## sigma - std deviation
## r - derivative index 
#
## Returns
## r-th derivative at x
###############################################################################

dnorm.deriv <- function(x, mu=0, sigma=1, deriv.order=0)
{
  r <- deriv.order
  phi <- dnorm(x, mean=mu, sd=sigma) 
  x <- (x - mu)
  
  if (r==0)
    return(phi)
  else if (r==1)
    derivt <- -x/sigma^2*phi
  else if (r==2)
    derivt <- (x^2-sigma^2)/sigma^4*phi
  else if (r==3)
    derivt <- -(x^3 - 3*x*sigma^2)/sigma^6*phi
  else if (r==4)
    derivt <- (x^4 - 6*x^2*sigma^2 + 3*sigma^4)/sigma^8*phi
  else if (r==5)
    derivt <- -(x^5 - 10*x^3*sigma^2 + 15*x*sigma^4)/sigma^10*phi
  else if (r==6)
    derivt <- (x^6 - 15*x^4*sigma^2 + 45*x^2*sigma^4 - 15*sigma^6)/sigma^12*phi
  else if (r==7)
    derivt <- -(x^7 - 21*x^5*sigma^2 + 105*x^3*sigma^4 - 105*x*sigma^6)/sigma^14*phi
  else if (r==8)
    derivt <- (x^8 - 28*x^6*sigma^2 + 210*x^4*sigma^4 - 420*x^2*sigma^6 + 105*sigma^8)/sigma^16*phi
  else if (r==9)
    derivt <- -(x^9 - 36*x^7*sigma^2 + 378*x^5*sigma^4 - 1260*x^3*sigma^6 + 945*x*sigma^8)/sigma^18*phi
  else if (r==10)
    derivt <- (x^10 - 45*x^8*sigma^2 + 630*x^6*sigma^4 - 3150*x^4*sigma^6 + 4725*x^2*sigma^8 - 945*sigma^10)/sigma^20*phi
  
  if (r > 10)
    stop ("Up to 10th order derivatives only")
    
  return(derivt)
}

###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation
#
## Parameters
## x - points to evaluate
## Sigma - variance matrix
## inc - 0 - exclude diagonals
##     - 1 - include diagonals
#
## Returns
## Double sum at x
###############################################################################

dnorm.sum <- function(x, sigma=1, inc=1, binned=FALSE, bin.par)
{
  d <- 1
  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, h=sigma)  
    n <- sum(bin.par$counts)
    fhatr <- kde.binned(bin.par=bin.par, h=sigma)
    sumval <- sum(bin.par$counts * n * fhatr$estimate)
    if (inc == 0) 
      sumval <- sumval - n*dnorm(x=0, mean=0, sd=sigma)
  }
  else
  {
    n <- length(x)
    if (n==1) sumval <- dnorm(x, mean=0, sd=sigma)
    else
    {  
      sumval <- 0
      for (i in 2:n)
      {sumval <- sumval + sum(dnorm(x=x[1:(i-1)]-x[i], mean=0, sd=sigma))}
      sumval<-2*sumval
      if (inc == 1) 
      sumval <- sumval + n*dnorm(x=0, mean=0, sd=sigma) 
    }
    
  }
 
  return(sumval)
}

dnorm.deriv.sum <- function(x, sigma, deriv.order, inc=1, binned=FALSE, bin.par, kfe=FALSE)
{
  r <- deriv.order
  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, h=sigma, supp=4+r)  

    fhatr <- kdde.binned(bin.par=bin.par, h=sigma, deriv.order=r)
    n <- sum(bin.par$counts)
    sumval <- sum(bin.par$counts * n * fhatr$estimate)
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, deriv.order=r)
  }
  else
  {
    n <- length(x)
    sumval <- 0
    for (i in 1:n)
      sumval <- sumval + sum(dnorm.deriv(x=x[i] - x, mu=0, sigma=sigma, deriv.order=r)) 
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, deriv.order=r)
  }

  
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))
  
  return(sumval)
  
}


###############################################################################
# Multivariate normal densities and derivatives
###############################################################################


###############################################################################
## Multivariate normal mixture - random sample
## 
## Parameters
## n - number of samples
## mus - matrix of means (each row is a vector of means from each component
##       density)
## Sigmas - matrix of covariance matrices (every d rows is a covariance matrix 
##          from each component density) 
## props - vector of mixing proportions 
## 
## Returns
## Vector of n observations from the normal mixture 
###############################################################################

rmvnorm.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), props=1, mixt.label=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")
  
  #if (is.vector(Sigmas))
  ##  return(rnorm.mixt(n=n, mus=mus, sigmas=Sigmas, props=props))
  
  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
   if (mixt.label)
     rand <- cbind(rmvnorm(n=n, mean=mus, sigma=Sigmas), rep(1, n))
   else
     rand <- cbind(rmvnorm(n=n, mean=mus, sigma=Sigmas))
    
  ### multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k)
    {
      ## compute random sample from normal mixture component
      if (n.prop[i] > 0)
      {       
        if (mixt.label)
          rand <- rbind(rand, cbind(rmvnorm(n=n.prop[i], mean=mus[i,], sigma=Sigmas[((i-1)*d+1) : (i*d),]), rep(i, n.prop[i])))
        else
          rand <- rbind(rand, rmvnorm(n=n.prop[i], mean=mus[i,], sigma=Sigmas[((i-1)*d+1) : (i*d),]))
      }    
    }
  }

  return(rand[sample(n),])
}


###############################################################################
## Multivariate normal mixture - density values
## 
## Parameters
## x - points to compute density at 
## mus - matrix of means
## Sigmas - matrix of covariance matrices 
## props - vector of mixing proportions 
## 
## Returns
## Density values from the normal mixture (at x)
###############################################################################

dmvnorm.mixt <- function(x, mus, Sigmas, props=1)
{  
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  if (is.vector(x)) d <- length(x)
  else d <- ncol(x)
  
  if (missing(mus)) mus <- rep(0,d)
  if (missing(Sigmas)) Sigmas <- diag(d)
  ##if (missing(deriv)) deriv <- rep(0,d)
   
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (is.matrix(mus)) mus <- mus[1,]
    dens <- dmvnorm(x=x, mean=mus, sigma=Sigmas[1:d,])
  }
  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0
    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dmvnorm(x, mean=mus[i,], sigma=Sigmas[((i-1)*d+1):(i*d),])
  }
  
  return(dens)
}   



###############################################################################
## Partial derivatives of the multivariate normal
## 
## Parameters
## x - points to evaluate at
## Sigma - variance
## r - derivative index 
#
## Returns
## r-th derivative at x
###############################################################################

dmvnorm.deriv.JEC <- function (x, mu, Sigma, deriv.order=0, Sdr.mat) #### The fastest one
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    n<-nrow(x)
    d<-ncol(x)
    if (missing(mu)) {
        mu <- rep(0, length = d)
    }
    if (missing(Sigma)) {
        Sigma <- diag(d)
    }
    r<-deriv.order
    if (missing(Sdr.mat)) {
        Sdr.mat <- Sdr(d=d,r=r)
    }
    
    ##### Normal density at x
    x.centred <- sweep(x, 2, mu)
    Sigmainv<- chol2inv(chol(Sigma))
    distval <- rowSums((x.centred %*% Sigmainv) * x.centred)  
    
    logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(d * log(2 * pi) + logdet + distval)/2
    dens<-matrix(exp(logretval),nrow=n)
    
    ##### Vector Hermite polynomial (Holmquist, 1996)
    
    if(r==0){ 
        result<-dens
        }

    if(r>0){
    vSigma<-vec(Sigma)
    Hr<-rep(0,d^r)
    ones<-rep(1,d^r)
    for(j in 0:floor(r/2)){
        cj<-(-1)^j*factorial(r)/(factorial(j)*2^j*factorial(r-2*j))
        vSigmaj<-Kpow(vSigma,j)
        xr2j<-mat.Kpow(A=x.centred,pow=r-2*j)
        Hr<-Hr+cj*mat.Kprod(U=xr2j,V=matrix(rep(vSigmaj,n),nrow=n,byrow=TRUE))
        }
        
    Sigmainvr<-Kpow(Sigmainv,r)        
    SirSdr<-Sigmainvr%*%Sdr.mat    
    Hr<-Hr%*%SirSdr
    result<-(-1)^r*(t(ones)%x%dens)*Hr
    }

    return(result)
}


dmvnorm.deriv <- function(x, mu, Sigma, deriv.order, Sdr.mat, duplicate=TRUE, add.index=FALSE, index.only=FALSE)
{
  r <- deriv.order
  sumr <- sum(r)
  
  if (missing(x)) d <- ncol(Sigma)
  else
  {
    if (is.vector(x))
      x <- t(as.matrix(x))
    
    if (is.data.frame(x)) x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
  }
  
  ## matrix of derivative indices
  
  if (sumr==0) ind.mat <- 0
  if (sumr==1) ind.mat <- diag(d)
  if (sumr==2) ind.mat <- K.sum(diag(d), diag(d))
  if (sumr==3) ind.mat <- K.sum(diag(d), K.sum(diag(d), diag(d)))
  if (sumr==4) ind.mat <- K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))
  if (sumr==5) ind.mat <- K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d)))))
  if (sumr==6) ind.mat <- K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))))
  if (sumr==7) ind.mat <- K.sum(diag(d), K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d)))))))
  if (sumr==8) ind.mat <- K.sum(diag(d), K.sum(diag(d), K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))))))

  if (index.only) return(ind.mat)

  ## derivatives
  
  if (missing(mu)) mu <- rep(0,d)
  if (missing(Sigma)) Sigma <- diag(d)
  if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d, r=sumr)
  ##for (i in 1:n) x[i,] <- x[i,] - mu
 
  mvh <- dmvnorm.deriv.JEC(x=x, mu=mu, Sigma=Sigma, deriv.order=sumr, Sdr.mat=Sdr.mat)

  if (length(r)>1)
  {
    mvh <- mvh[,!duplicated(ind.mat)]
    deriv.ind <- unique(ind.mat)
    which.deriv <- which.mat(r, deriv.ind)
    if (is.vector(mvh)) return(mvh[which.deriv])
    else return(mvh[,which.deriv])
  }
  else
  {
    if (r>1)
    {  
      if (!duplicate)
      { 
        mvh <- mvh[,!duplicated(ind.mat)]
        ind.mat <- unique(ind.mat)
      }
    }
    ##else
    ##  mvh <- (-1)^sumr*mvh

    if (add.index) return(list(deriv=mvh, deriv.ind=ind.mat))
    else return(deriv=mvh)
  }
}




###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation
#
## Parameters
## x - points to evaluate
## Sigma - variance matrix
## inc - 0 - exclude diagonals
##     - 1 - include diagonals
#
## Returns
## Double sum at x
##############################################################################

dmvnorm.sum <- function(x, Sigma, inc=1, binned=FALSE, bin.par, diff=FALSE)
{
  if (binned)
  {
    if (!is.diagonal(Sigma))
      stop("Binned estimation defined for diagonal Sigma only")
    if (missing(bin.par)) bin.par <- binning(x, H=Sigma)  
    
    fhatr <- kde.binned(bin.par=bin.par, H=Sigma)$estimate 
    n <- sum(bin.par$counts)
    d <- ncol(Sigma)
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm(x=rep(0,d), mean=rep(0,d), sigma=Sigma)
  }
  else
  {
    d <- ncol(Sigma)

    if (d>=2)
    {
      if(!diff)
      {
        n <- nrow(x)
        d <- ncol(x)
        difs <- differences(x, upper=TRUE)
      }
      else
      {
        n <- (-1 + sqrt(1+8*nrow(x)))/2
        d <- ncol(x)
        difs <- x
      }
      sumval <- sum(dmvnorm(difs, mean=rep(0,d), sigma=Sigma))
      
      if (inc==0)
        sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), mean=rep(0,d), sigma=Sigma)
      if (inc==1)
        sumval <- 2*sumval - n*dmvnorm(rep(0,d), mean=rep(0,d), sigma=Sigma) 
    }
  
  }
  return(sumval)
}


dmvnorm.deriv.sum <- function(x, Sigma, deriv.order=0, inc=1, binned=FALSE, bin.par, diff=FALSE, kfe=FALSE)
{
  r <- deriv.order
  if (binned)
  {
    if (!is.diagonal(Sigma))
        stop("Binned estimation defined for diagonal Sigma only")
    if (missing(bin.par)) bin.par <- binning(x, H=Sigma)  
    
    fhatr <- kdde.binned(bin.par=bin.par, H=Sigma, deriv.order=r)$estimate 
    n <- sum(bin.par$counts)
    d <- ncol(Sigma)
    sumval <- sum(bin.par$counts * n * fhatr)
  }
  else
  {
    if(!diff)
    {
      ## can still optimise this for even order derivatives
      n <- nrow(x)
      d <- ncol(x)
      difs <- differences(x, upper=FALSE)
    }
    else
    {
      n <- sqrt(nrow(x)) #(-1 + sqrt(1+8*nrow(x)))/2
      d <- ncol(x)
      difs <- x
    }

    mvh.temp <- dmvnorm.deriv(difs, mu=rep(0,d), Sigma=Sigma, deriv.order=r)
    if (is.matrix(mvh.temp))
      sumval <- apply(mvh.temp,2, sum)
    else
      sumval <- sum(mvh.temp)
  }
  
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=Sigma, deriv=r)
  
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))
  
  return(sumval)
}


## Used in Hbcv

dmvnorm.deriv.2d.xxt.sum <- function(x, Sigma, r)
{
  if (is.vector(x))
  {
    n <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnormd4_2d_xxt_sum", as.double(x1), as.double(x2),
               as.double(viSigma), as.integer(r), as.integer(n), 
               as.double(rep(0,4)), PACKAGE="ks")
  
  ## above C functions dmvnorm4_2d_xxt_sum only computes the upper triangular
  ## half so need to reflect along the diagonal to compute whole sum
  sumval <- 2*result[[6]]
  
  return(invvec(sumval))
} 


###############################################################################
## Compute moments of multivariate normal mixture
###############################################################################

moments.mixt <- function (mus, Sigmas, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")
  d <- ncol(Sigmas)
  k <- length(props)
  mn <- rep(0, d)
  va <- matrix(0, nr=d, nc=d)
  for (i in 1:k)
  {
    mn <- mn + props[i] * mus[i,]
    va <- va + props[i] * (Sigmas[((i-1)*d+1):(i*d),] + mus[i,] %*% t(mus[i,]))
  } 
  va <- va + mn %*% t(mn)
  return( list(mean=mn, var=va))
}


###############################################################################
## Creates plots of mixture density functions
#
## Parameters
## mus - means
## Sigmas - variances
## props - vector of proportions of each mixture component 
## dfs - degrees of freedom
## dist - "normal" - normal mixture
##      - "t" - t mixture
## ...
###############################################################################


plotmixt <- function(mus, Sigmas, props, dfs, dist="normal", ...)
{
  if (ncol(Sigmas)==2)
    plotmixt.2d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, ...)
  else if (ncol(Sigmas)==3)
    plotmixt.3d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, ...) 
}


plotmixt.2d <- function(mus, Sigmas, props, dfs, dist="normal",
    xlim, ylim, gridsize, display="slice", cont=c(25,50,75), abs.cont,
    lty, xlab="x", ylab="y", zlab="Density function",
    theta=-30, phi=40, d=4, add=FALSE, drawlabels=TRUE, nrand=1e5, ...)
{
  dist <- tolower(substr(dist,1,1))
  maxSigmas <- 4*max(Sigmas)

  if (is.vector(mus))
    mus <- as.matrix(t(mus))

  if (missing(xlim))
    xlim <- c(min(mus[,1]) - maxSigmas, max(mus[,1]) + maxSigmas)
  if (missing(ylim))
    ylim <- c(min(mus[,2]) - maxSigmas, max(mus[,2]) + maxSigmas)

  if (missing(gridsize))
    gridsize <- rep(51,2)
              
  x <- seq(xlim[1], xlim[2], length=gridsize[1])
  y <- seq(ylim[1], ylim[2], length=gridsize[2])
  xy <- permute(list(x, y))

  d <- ncol(Sigmas)
  
  if (dist=="n")
    dens <- dmvnorm.mixt(xy, mu=mus, Sigma=Sigmas, props=props)
  
  else if (dist=="t")
    dens <- dmvt.mixt(xy, mu=mus, Sigma=Sigmas, props=props, dfs=dfs)

  dens.mat <- matrix(dens, nc=length(x), byrow=FALSE)
   
  disp <- substr(display,1,1)

  if (disp=="p")
    persp(x, y, dens.mat, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab,
          zlab=zlab, ...)

  else if (disp=="s")
  {
    if (dist=="n")
    {
      x.rand <- rmvnorm.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props)
      dens.rand <- dmvnorm.mixt(x.rand, mus=mus, Sigmas=Sigmas, props=props)
    }
    else if (dist=="t")
    {
      x.rand <- rmvt.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
      dens.rand <- dmvt.mixt(x.rand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
    }
    
    if (missing(lty))
      lty <- 1
    if (missing(abs.cont))
      hts <- quantile(dens.rand, prob=(100 - cont)/100)
    else
      hts <- abs.cont
    
    if (!add)
      plot(x, y, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
    
    
    for (i in 1:length(hts)) 
    {
      scale <- cont[i]/hts[i]
      if (missing(abs.cont))
        contour(x, y, dens.mat*scale, level=hts[i]*scale, add=TRUE,
                drawlabels=drawlabels, lty=lty, ...)
      else
        contour(x, y, dens.mat, level=hts[i], add=TRUE,
                drawlabels=drawlabels, lty=lty, ...)
    }
  }
  
  else if (disp=="i")
    image(x, y, dens.mat, xlab=xlab, ylab=ylab, ...)
  else if (disp=="f")
    filled.contour(x, y, dens.mat, xlab=xlab, ylab=ylab, ...)
    
}

plotmixt.3d <- function(mus, Sigmas, props, dfs, cont=c(25,50,75), abs.cont,
    dist="normal", xlim, ylim, zlim, gridsize, alphavec, colors, add=FALSE, nrand=1e5, ...)
{
  require(rgl)
  require(misc3d)
  d <- 3
  dist <- tolower(substr(dist,1,1))
  maxsd <- sqrt(apply(Sigmas, 2, max))

  if (is.vector(mus))
    mus <- as.matrix(t(mus))

  if (missing(xlim))
    xlim <- c(min(mus[,1]) - 4*maxsd[1], max(mus[,1]) + 4*maxsd[1])
  if (missing(ylim))
    ylim <- c(min(mus[,2]) - 4*maxsd[2], max(mus[,2]) + 4*maxsd[2])
  if (missing(zlim))
    zlim <- c(min(mus[,3]) - 4*maxsd[3], max(mus[,3]) + 4*maxsd[3])

  if (missing(gridsize))
    gridsize <- rep(51,d)
  
  x <- seq(xlim[1], xlim[2], length=gridsize[1])
  y <- seq(ylim[1], ylim[2], length=gridsize[2])
  z <- seq(zlim[1], zlim[2], length=gridsize[3])
  xy <- permute(list(x,y))

  dens.array <- array(0, dim=gridsize)
  
  for (i in 1:length(z))
  {
    if (dist=="n")
      dens <- dmvnorm.mixt(cbind(xy, z[i]), mu=mus, Sigma=Sigmas, props=props)
    else if (dist=="t")
      dens <- dmvt.mixt(cbind(xy, z[i]), mu=mus, Sigma=Sigmas, dfs=dfs, props=props)
    
    dens.mat <- matrix(dens, nc=length(x), byrow=FALSE)
    dens.array[,,i] <- dens.mat
  }
  
  if (dist=="n")
  {  
    x.rand <- rmvnorm.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props)
    dens.rand <- dmvnorm.mixt(x.rand, mus=mus, Sigmas=Sigmas, props=props)
  }
  else if (dist=="t")
  {
    x.rand <- rmvt.mixt(n=nrand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
    dens.rand <- dmvt.mixt(x.rand, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)
  }
  if (missing(abs.cont))
    hts <- quantile(dens.rand, prob = (100 - cont)/100)
  else
    hts <- abs.cont

  nc <- length(hts)
  if (missing(colors))
    colors <- rev(heat.colors(nc))
  
  if (missing(alphavec))
    alphavec <- seq(0.1,0.5,length=nc)

  ##plot3d(x, y, z, type="n", add=add, ...)
  if (!add) clear3d()
  
  for (i in 1:nc) 
  {
    contour3d(dens.array, level=hts[nc-i+1], x,y,z, add=TRUE, color=colors[i],
             alpha=alphavec[i],...)
  }
  decorate3d(...)
  fhat <- list(eval.points=list(x, y, z), estimate=dens.array, cont=hts)
 
  invisible(fhat)
}


###############################################################################
## Multivariate t - density values
#
## Parameters
## x - points to compute density     
## mu - vector of means 
## Sigma - dispersion matrix
## df - degrees of freedom
#
## Returns
## Value of multivariate t density at x
###############################################################################

dmvt <- function(x, mu, Sigma, df)
{   
  if(is.vector(x))
    x <- matrix(x, ncol=length(x))
  d <- ncol(Sigma)
  detSigma <- det(Sigma) 
  dens <- (1+ mahalanobis(x, center=mu, cov=Sigma)/df)^(-(d+df)/2)
  dens <- dens * gamma((df+d)/2) / ((df*pi)^(d/2)*gamma(df/2)*detSigma^(1/2))
  
  return(dens)
}


###############################################################################
## Multivariate t mixture - density values
#
## Parameters
## x - points to compute density at    
## mus - vector of means 
## Sigmas - dispersion matrices
## dfs - degrees of freedom
## props - vector of mixing proportions
#
## Returns
## Value of multivariate t mixture density at x
###############################################################################

dmvt.mixt <- function(x, mus, Sigmas, dfs, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")
  
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dmvt(x, mu=mus, Sigma=Sigmas, df=dfs)
  
  ## multiple component mixture
  else   
  {   
    if (is.vector(mus)) d <- length(mus)
    else d <- ncol(mus)
    k <- length(props)
    dens <- 0      
    for (i in 1:k)
      dens <- dens+props[i]*dmvt(x,mu=mus[i,],Sigma=Sigmas[((i-1)*d+1):(i*d),],
                                 df=dfs[i])
  }
  
  return(dens)
}


###############################################################################
## Multivariate t mixture - random sample
## 
## Parameters
## n - number of samples
## mus - means 
## Sigmas - matrix of dispersion matrices
## dfs - vector of degrees of freedom
## props - vector of mixing proportions 
## 
## Returns
## Vector of n observations from the t mixture
###############################################################################

rmvt.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))  
    stop("Proportions don't sum to one\n")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")  

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    rand <- rmvt(n=n, sigma=Sigmas, df=dfs)
    for (i in 1:length(mus))
      rand[,i] <- rand[,i] + mus[i]
  }
  
  ## multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    ## compute number to be drawn from each component 
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))

    ## generate random samples from each component
    rand <- numeric(0)  
    for (i in 1:k)
    {
      if (n.prop[i] > 0)
      {  
        rand.temp<-rmvt(n=n.prop[i],sigma=Sigmas[((i-1)*d+1):(i*d),],df=dfs[i])
        for (j in 1:length(mus[k,]))
          rand.temp[,j] <- rand.temp[,j] + mus[i,j]
       
        rand <- rbind(rand, rand.temp)
      }
    }
  }
  
  return(rand[sample(n),])
}


