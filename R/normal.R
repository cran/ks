
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
    dens <- dnorm(x, mean=mus, sd=sigmas)

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

dnorm.deriv <- function(x, mu=0, sigma=1, r=0)
{
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
    ##sumval <- bkfe(x=bin.par$counts, bandwidth=sigma, drv=0, binned=TRUE, range.x=range(bin.par$eval.points))
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, r=0, sigma=sigma)
  }
  else
  {
    n <- length(x)
    sumval <- 0
    for (i in 1:n)
      sumval <- sumval + sum(dnorm(x=x[i] - x, mean=0, sd=sigma))
    
    if (inc == 0) 
      sumval <- sumval - n*dnorm(x=0, mean=0, sd=sigma)
   } 
  
  return(sumval)
}

dnorm.deriv.sum <- function(x, sigma, r, inc=1, binned=FALSE, bin.par, kfe=FALSE)
{
  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, h=sigma, supp=4+r)  

    fhatr <- kdde.binned(bin.par=bin.par, h=sigma, r=r)
    n <- sum(bin.par$counts)
    sumval <- sum(bin.par$counts * n * fhatr$estimate)
    ##sumval <- n*bkfe(x=bin.par$counts, bandwidth=sigma, drv=r, binned=TRUE, range.x=range(bin.par$eval.points))
  }
  else
  {
    n <- length(x)
    sumval <- 0
    for (i in 1:n)
      sumval <- sumval + sum(dnorm.deriv(x=x[i] - x, mu=0, sigma=sigma, r=r)) 
  }

  if (inc == 0) 
    sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, r=r)

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
    dens <- dmvnorm(x=x, mean=mus, sigma=Sigmas)
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

### for diagonal Sigma
dmvnorm.deriv.diag <- function(x, mu, Sigma, r)
{
  if (is.vector(x))
    x <- t(as.matrix(x))

  if (is.data.frame(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  
  if (missing(mu)) mu <- rep(0,d)
  if (missing(Sigma)) Sigma <- diag(d)

  for (i in 1:n)
    x[i,] <- x[i,] - mu  
  
  if (d==2) return(dmvnorm.deriv.2d(x=x, Sigma=Sigma, r=r))
  if (d==3) return(dmvnorm.deriv.3d(x=x, Sigma=Sigma, r=r))
  if (d==4) return(dmvnorm.deriv.4d(x=x, Sigma=Sigma, r=r))
  if (d==5) return(dmvnorm.deriv.5d(x=x, Sigma=Sigma, r=r))
  if (d==6) return(dmvnorm.deriv.6d(x=x, Sigma=Sigma, r=r))

}

### for general Sigma
dmvnorm.deriv <- function(x, mu, Sigma, r, Sdr.mat)
{   
  if (is.vector(x))
    x <- t(as.matrix(x))

  if (is.data.frame(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  
  sumr <- sum(r)
  
  if (missing(mu)) mu <- rep(0,d)
  if (missing(Sigma)) Sigma <- diag(d)

  for (i in 1:n)
      x[i,] <- x[i,] - mu
  
  if (missing(Sdr.mat) & d >=2)
    Sdr.mat <- Sdr(d=d, r=sumr)

  dens <- dmvnorm(x=x, mean=mu, sigma=Sigma) 
  vSigma <- vec(Sigma)
  Sigmainv <- chol2inv(chol(Sigma))
  mvh <- matrix(0, nrow=n, ncol=d^sumr)
 
  if (sumr==0)
    mvh <- dens
  if (sumr==1)
    mvh <- x*dens
  if (sumr==2)
  {
    Sinv <- (Sigmainv %x% Sigmainv) %*%  Sdr.mat
    for (i in 1:n)
      mvh[i,]  <-  Sinv %*% ((x[i,] %x% x[i,]) - vSigma) * dens[i]
    ind.mat <- K.sum(diag(d), diag(d))
  }

  if (sumr==3)
  {
    Sinv <- K.pow(Sigmainv,3) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,], 3) - 3*x[i,] %x% vSigma) * dens[i]
    ind.mat <- K.sum(diag(d), K.sum(diag(d), diag(d)))
  }

  if (sumr==4)
  {    
    Sinv <- K.pow(Sigmainv,4) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,], 4) - 6*K.pow(x[i,],2) %x% vSigma + 3*vSigma %x% vSigma) * dens[i]
    ind.mat <- K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))
  }

  if (sumr==5)
  {
    Sinv <- K.pow(Sigmainv,5) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,], 5) - 10*K.pow(x[i,],3) %x% vSigma + 15*x[i,] %x% K.pow(vSigma,2)) * dens[i]
    ind.mat <- K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d)))))
  }

  if (sumr==6)
  {
    Sinv <- K.pow(Sigmainv,6) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,],6) - 15*K.pow(x[i,],4) %x% vSigma + 45*K.pow(x[i,],2) %x% K.pow(vSigma,2) - 15*K.pow(vSigma,3)) * dens[i]
    ind.mat <- K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))))
  }

  if (sumr==7)
  {
    Sinv <- K.pow(Sigmainv,7) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,],7) - 21*K.pow(x[i,],5) %x% vSigma + 105*K.pow(x[i,],3) %x% K.pow(vSigma,2)) * dens[i]
    ind.mat <- K.sum(diag(d), K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d)))))))
  }

  if (sumr==8)
  {
    Sinv <- K.pow(Sigmainv,8) %*% Sdr.mat
    for (i in 1:n)
      mvh[i,] <- Sinv %*% (K.pow(x[i,],8) - 28*K.pow(x[i,],6) %x% vSigma + 210*K.pow(x[i,],4) %x% K.pow(vSigma,2) - 420*K.pow(x[i,],2) %x% K.pow(vSigma,3) + 105*K.pow(vSigma, 4)) * dens[i]
    ind.mat <- K.sum(diag(d), K.sum(diag(d), K.sum(diag(d),K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), K.sum(diag(d), diag(d))))))))
  }

  if (sumr > 8)
    stop ("Up to 8th order derivatives only")

  mvh <- (-1)^sumr*mvh[,!duplicated(ind.mat)]

  deriv.ind <- unique(ind.mat)
  if (length(r)>1)
  {
    which.deriv <- which.mat(r, deriv.ind)
    if (is.vector(mvh)) return(mvh[which.deriv])
    else return(mvh[,which.deriv])
  }
  else
    return(list(deriv=mvh, deriv.ind=deriv.ind))
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
    if (!identical(diag(diag(Sigma)), Sigma))
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
    ### Need to rewrite this for d <=6
    d <- ncol(Sigma)
    if (d==2) sumval <- dmvnorm.2d.sum(x=x, Sigma=Sigma, inc=inc)
    if (d==3) sumval <- dmvnorm.3d.sum(x=x, Sigma=Sigma, inc=inc)
    if (d==4) sumval <- dmvnorm.4d.sum(x=x, Sigma=Sigma, inc=inc)
    if (d==5) sumval <- dmvnorm.5d.sum(x=x, Sigma=Sigma, inc=inc)
    if (d==6) sumval <- dmvnorm.6d.sum(x=x, Sigma=Sigma, inc=inc)
    
    if (d>6)
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


dmvnorm.deriv.sum <- function(x, Sigma, r, inc=1, binned=FALSE, bin.par, diff=FALSE, kfe=FALSE)
{
  if (binned)
  {
    if (!is.diagonal(Sigma))
      stop("Binned estimation defined for diagonal Sigma only")
    if (missing(bin.par)) bin.par <- binning(x, H=Sigma)  

    fhatr <- kdde.binned(bin.par=bin.par, H=Sigma, r=r)$estimate 
    n <- sum(bin.par$counts)
    d <- ncol(Sigma)
    sumval <- sum(bin.par$counts * n * fhatr)
  }
  else
  {
    if(!diff)
    {
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
    if (is.diagonal(Sigma))
      sumval <- sum(dmvnorm.deriv.diag(difs, mu=rep(0,d), Sigma=Sigma, r=r))
    else
      sumval <- sum(dmvnorm.deriv(difs, mu=rep(0,d), Sigma=Sigma, r=r))
  }

  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=Sigma, r=r)
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))
  
  return(sumval)
}









###############################################################################
## Partial derivatives of the bivariate normal (mean 0) 
## 
## Parameters
## x - points to evaluate at
## Sigma - variance matrix
## r - (r1, r2) vector of partial derivative indices 
#
## Returns
## r-th partial derivative at x
###############################################################################

dmvnorm.deriv.2d <- function(x, Sigma, r) 
{
  if (is.vector(x))
  {    
    n <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else 
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }

  if (sum(r) == 0)
    return(dmvnorm(x, c(0,0), Sigma))
  ## first order derivatives  
  else if (sum(r) == 1)
    derivt <- .C("dmvnormd1_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  ## third order derivatives
  else if (sum(r) == 2)
    derivt <- .C("dmvnormd2_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  ## second order derivatives
  else if (sum(r) == 3)
    derivt <- .C("dmvnormd3_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  ## fourth order derivatives
  else if (sum(r) == 4)
    derivt <- .C("dmvnormd4_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  ## fifth order derivatives
  else if (sum(r) == 5)
    derivt <- .C("dmvnormd5_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  ## sixth order derivatives
  else if (sum(r) == 6)
    derivt <- .C("dmvnormd6_2d", as.double(x1), as.double(x2), 
            as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  else
    stop("Only works for up to 6th order partial derivatives")
  
  return(derivt[[6]])
}         

###############################################################################
## Partial derivatives of the 3-variate normal (mean 0, diag variance matrix) 
##############################################################################

dmvnorm.deriv.3d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment
  
  ##d <- 3
  if (sum(r) > 6)
    stop("Only works for up to 6th order partial derivatives")

  if (is.vector(x))
  {    
    ##n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; 
  }
  else 
  {
    ##n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; 
  }
  
  y1 <- cbind(x1,x2)
  y2 <- x3
  r1 <- r[c(1,2)]
  r2 <- r[3]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- Sigma[3,3]

  ## Use existing 2-dim derivatives to compute 3-dim derivative for diag
  ## variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dnorm.deriv(x=y2, mu=0, sigma=sqrt(Sigma2), r=r2)
  derivt <- derivt1*derivt2

  return(derivt)
}         

##############################################################################
## Partial derivatives of the 4-variate normal (mean 0, diag variance matrix) 
##############################################################################

dmvnorm.deriv.4d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment
  
  ##d <- 4
  if (sum(r) > 6)
    stop("Only works for up to 6th order partial derivatives")

  if (is.vector(x))
  {    
    ##n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
  }
  else 
  {
    ##n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
  }
  
  y1 <- cbind(x1,x2)
  y2 <- cbind(x3,x4)
  r1 <- r[c(1,2)]
  r2 <- r[c(3,4)]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- diag(c(Sigma[3,3], Sigma[4,4]))

  ## Use existing 2-dim derivatives to compute 4-dim derivative for diag
  ## variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt <- derivt1*derivt2

  return(derivt)
}
###############################################################################
## Partial derivatives of the 5-variate normal (mean 0, diag variance matrix) 
## 
## Parameters
## x - points to evaluate at
## Sigma - variance matrix
## r - vector of partial derivative indices 
#
## Returns
## r-th partial derivative at x
##############################################################################

dmvnorm.deriv.5d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment

  ##d <- 5
  if (sum(r) > 6)
    stop("Only works for 2nd, 4th and 6th order partial derivatives")
 
  if (is.vector(x))
  {    
    ##n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]; 
  }
  else 
  {
    ##n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4]; x5 <- x[,5]; 
  }
  
  y1 <- cbind(x1,x2)
  y2 <- cbind(x3,x4)
  y3 <- x5
  r1 <- r[c(1,2)]
  r2 <- r[c(3,4)]
  r3 <- r[5]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- diag(c(Sigma[3,3], Sigma[4,4]))
  Sigma3 <- Sigma[5,5]
  
  ## Use existing 2-dim derivatives to compute 5-dim derivative for diag
  ## variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt3 <- dnorm.deriv(x=y3, mu=0, sigma=sqrt(Sigma3), r=r3)
  derivt <- derivt1*derivt2*derivt3

  return(derivt)
}         


###############################################################################
## Partial derivatives of the 6-variate normal (mean 0, diag variance matrix) 
## 
## Parameters
## x - points to evaluate at
## Sigma - variance matrix
## r - vector of partial derivative indices 
#
## Returns
## r-th partial derivative at x
##############################################################################

dmvnorm.deriv.6d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment

  ##d <- 6
  if (sum(r) > 6)
    stop("Only works for 2nd, 4th and 6th order partial derivatives")

  if (is.vector(x))
  {    
    ##n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]; x6 <- x[6];
  }
  else 
  {
    ##n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4]; x5 <- x[,5]; x6 <- x[,6];
  }
  
  y1 <- cbind(x1,x2)
  y2 <- cbind(x3,x4)
  y3 <- cbind(x5,x6)
  r1 <- r[c(1,2)]
  r2 <- r[c(3,4)]
  r3 <- r[c(5,6)]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- diag(c(Sigma[3,3], Sigma[4,4]))
  Sigma3 <- diag(c(Sigma[5,5], Sigma[6,6]))
  
  ## Use existing 2-dim derivatives to compute 6-dim derivative for diag
  ## variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt3 <- dmvnorm.deriv.2d(y3, Sigma3, r3)
  derivt <- derivt1*derivt2*derivt3

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

dmvnorm.2d.sum <- function(x, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); d <- 2; x1 <- x[,1]; x2 <- x[,2]
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_2d_sum", as.double(x1), as.double(x2),
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[6]]
  
  ## above C function mvnorm_2d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(c(0,0), mean=c(0,0), sigma=Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(c(0,0), mean=c(0,0), sigma=Sigma) 
  
  return(sumval)
}

dmvnorm.deriv.2d.sum <- function(x, Sigma, r, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
     
  if (sum(r)==0)
    return(dmvnorm.2d.sum(x, Sigma, inc=inc))
  ## 1st order derivatives
  else if (sum(r)==1)
    derivt <- .C("dmvnormd1_2d_sum", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0), PACKAGE="ks")
  ## 2nd order derivativeselse if (sum(r)==2)
  else if (sum(r)==2)
    derivt <- .C("dmvnormd2_2d_sum", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0), PACKAGE="ks")
  ## 3rd order derivativeselse if (sum(r)==3)
  else if (sum(r)==3)
    derivt <- .C("dmvnormd3_2d_sum", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0), PACKAGE="ks")
  ## fourth order derivatives
  else if (sum(r) == 4)
    derivt <- .C("dmvnormd4_2d_sum", as.double(x1), as.double(x2), 
                   as.double(vec(Sigma)), as.integer(r), as.integer(n),
                 as.double(0), PACKAGE="ks")
  ## fifth order derivatives
  else if (sum(r) == 5)
    derivt <- .C("dmvnormd5_2d_sum", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n),
                 as.double(0), PACKAGE="ks")
  ## sixth order derivatives
  else if (sum(r) == 6)
    derivt <- .C("dmvnormd6_2d_sum", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n),
                   as.double(0), PACKAGE="ks")
  else
    stop("Only works for up to 6th order partial derivatives")
  
  if (sum(r)==0)
    sumval <- derivt
    else
      sumval <- derivt[[6]]
  
  ## above C functions mvnorm?_2d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma=Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma=Sigma) 

  return(sumval)
}         


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
## Double sum  of K(X_i - X_j) used in density derivative estimation - 3-dim
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

dmvnorm.3d.sum <- function(x, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3];
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_3d_sum", as.double(x1),as.double(x2), as.double(x3), 
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[7]]
  
  ## above C function mvnorm_3d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), sigma=Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), sigma=Sigma) 
  
  return(sumval)
}

dmvnorm.deriv.3d.sum <- function(x, Sigma, r, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; x1 <- x[1]; x2 <- x[2] ; x3 <- x[3];  
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3];   
  }

  d <- 3
  sumval <- 0

  if (binned)
  {
    ## drvkde computes include-diagonals estimate
    fhatr <- drvkde(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.3d(x=rep(0,d), r=r, Sigma=Sigma)
  }
  else
  {  
    for (j in 1:n)
    {  
      y1 <- x1 - x1[j]
      y2 <- x2 - x2[j]
      y3 <- x3 - x3[j]
      sumval <- sumval + sum(dmvnorm.deriv.3d(cbind(y1, y2, y3), Sigma, r))
    }
    
    if (inc==0)
      sumval <- sumval - n*dmvnorm.deriv.3d(rep(0,d), Sigma, r)
  }

  
  return(sumval)
}

###############################################################################
# Double sum  of K(X_i - X_j) used in density derivative estimation - 4-dim
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

dmvnorm.4d.sum <- function(x, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
  }

  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_4d_sum", as.double(x1), as.double(x2),
               as.double(x3), as.double(x4),
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[8]]
    
  ## above C function mvnorm_2d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), sigma=Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), sigma=Sigma) 

  return(sumval)
}

dmvnorm.deriv.4d.sum <- function(x, Sigma, r, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; x1 <- x[1]; x2 <- x[2] ;x3 <- x[3]; x4 <- x[4]; 
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
  }

  d <- 4
  sumval <- 0

  if (binned)
  {
    ## drvkde computes include-diagonals estimate
    fhatr <- drvkde(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.4d(x=rep(0,d), r=r, Sigma=Sigma)
  }
  else
  {  
    for (j in 1:n)
    {  
      y1 <- x1 - x1[j]
      y2 <- x2 - x2[j]
      y3 <- x3 - x3[j]
      y4 <- x4 - x4[j]
      sumval <- sumval + sum(dmvnorm.deriv.4d(cbind(y1, y2, y3, y4), Sigma=Sigma, r=r))
    }
  
    if (inc==0)
      sumval <- sumval - n*dmvnorm.deriv.4d(rep(0,d), Sigma, r)
  }
  
  return(sumval)
}


###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation - 5-dim
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

dmvnorm.5d.sum <- function(x, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
    x5 <- x[5]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
    x5 <- x[,5]; 
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_5d_sum", as.double(x1), as.double(x2),
               as.double(x3), as.double(x4), as.double(x5),
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[9]]

  ## above C function mvnorm_5d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  
  return(sumval)
}

dmvnorm.deriv.5d.sum <- function(x, Sigma, r, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
    x5 <- x[5]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
    x5 <- x[,5]; 
  }

  d <- 5
  sumval <- 0

  for (j in 1:n)
  {  
    y1 <- x1 - x1[j]
    y2 <- x2 - x2[j]
    y3 <- x3 - x3[j]
    y4 <- x4 - x4[j]
    y5 <- x5 - x5[j]
    sumval <- sumval + sum(dmvnorm.deriv.5d(cbind(y1, y2, y3, y4, y5),Sigma,r))
  }
  
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv.5d(rep(0,d), Sigma, r)
    
  return(sumval)
}

###############################################################################
## Double sum  of K(X_i - X_j) used in density derivative estimation - 6-dim
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

dmvnorm.6d.sum <- function(x, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
    x5 <- x[5]; x6 <- x[6]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
    x5 <- x[,5]; x6 <- x[,6];  
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_6d_sum", as.double(x1), as.double(x2),
               as.double(x3), as.double(x4), as.double(x5), as.double(x6),
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[10]]

  ## above C function mvnorm_6d_sum only computes the upper triangular half
  ## so need to reflect along the diagonal and then subtract appropriate
  ## amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  
  return(sumval)
}

dmvnorm.deriv.6d.sum <- function(x, Sigma, r, inc=1)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
    x5 <- x[5]; x6 <- x[6]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
    x5 <- x[,5]; x6 <- x[,6];
  }

  d <- 6
  sumval <- 0

  for (j in 1:n)
  {  
    y1 <- x1 - x1[j]
    y2 <- x2 - x2[j]
    y3 <- x3 - x3[j]
    y4 <- x4 - x4[j]
    y5 <- x5 - x5[j]
    y6 <- x6 - x6[j]
    sumval <- sumval + sum(dmvnorm.deriv.6d(cbind(y1, y2, y3, y4, y5, y6),Sigma,r))
  }
  
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv.6d(rep(0,d), Sigma, r)
    
  return(sumval)
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
  maxSigmas <- 3.7*max(Sigmas)

  if (is.vector(mus))
    mus <- as.matrix(t(mus))
  
  if (missing(xlim))
    xlim <- c(min(mus[,1]) - maxSigmas, max(mus[,1]) + maxSigmas)
  if (missing(ylim))
    ylim <- c(min(mus[,2]) - maxSigmas, max(mus[,2]) + maxSigmas)
  if (missing(zlim))
    zlim <- c(min(mus[,3]) - maxSigmas, max(mus[,3]) + maxSigmas)
  
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

  plot3d(x, y, z, type="n", add=add, ...)

  for (i in 1:nc) 
  {
    ##scale <- cont[i]/hts[i]
    contour3d(dens.array, level=hts[nc-i+1],x, y, z, add=TRUE, color=colors[i],
             alpha=alphavec[i], ...)
  }
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
        rand.temp<-rmvt(n=n.prop[i],sigma=Sigmas[((i-1)*d+1):(i*d),],df=dfs[k])
        for (j in 1:length(mus[k,]))
          rand.temp[,j] <- rand.temp[,j] + mus[i,j]
       
        rand <- rbind(rand, rand.temp)
      }
    }
  }
  
  return(rand[sample(n),])
}


