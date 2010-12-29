
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


dnorm.deriv.sum <- function(x, sigma, deriv.order, inc=1, binned=FALSE, bin.par, kfe=FALSE)
{
  r <- deriv.order
  n <- length(x)
  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, h=sigma, supp=4+r) 
    est <- kdde.binned(x=x, H=sigma^2, h=sigma, deriv.order=r, bin.par=bin.par)$estimate
    sumval <- sum(bin.par$counts*est*n)
    if (inc == 0) 
      sumval <- sumval - n*dnorm.deriv(x=0, mu=0, sigma=sigma, deriv.order=r)
  }
  else
  {
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

dnorm.deriv.mixt <- function(x, mus=0, sigmas=1, props=1, deriv.order=0)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm.deriv(x, mu=mus[1], sigma=sigmas[1], deriv.order=deriv.order)

  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0

    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dnorm.deriv(x=x, mu=mus[i], sigma=sigmas[i], deriv.order=deriv.order)
  }
  
  return(dens)
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


dmvnorm.deriv <- function(x, mu, Sigma, deriv.order, Sdr.mat, deriv.vec=TRUE, add.index=FALSE, only.index=FALSE)
{
  r <- deriv.order
  if(length(r)>1) stop("deriv.order should be a non-negative integer.")
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
  ind.mat <- 0
  sumr.counter <- sumr
  if (sumr>=1) ind.mat <- diag(d)
  {
    while (sumr.counter >1)
    {
      ind.mat <- Ksum(diag(d), ind.mat)
      sumr.counter <- sumr.counter - 1
    }
  }
  ind.mat.minimal <- unique(ind.mat)
  ind.mat.minimal.logical <- !duplicated(ind.mat)
  
  if (only.index)
    if (deriv.vec) return (ind.mat)
    else return(ind.mat.minimal)
  
  ## compute derivatives

  if (missing(mu)) mu <- rep(0,d)
  if (missing(Sigma)) Sigma <- diag(d)
  if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d, r=sumr)

  ## Code by Jose Chacon 
  ## Normal density at x
  
  x.centred <- sweep(x, 2, mu)
  Sigmainv<- chol2inv(chol(Sigma))
  distval <- rowSums((x.centred %*% Sigmainv) * x.centred)  
  
  logdet <- sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(d * log(2 * pi) + logdet + distval)/2
  
  ## Vector Hermite polynomial (Holmquist, 1996)
    
  if(r==0){ 
    mvh <- matrix(exp(logretval), nrow=n)
  }

  
  if(r>0){
    vSigma<-vec(Sigma)
    ones <- rep(1,d^r)
    Sigmainvr <- Kpow(Sigmainv,r)        
    SirSdr <- Sigmainvr %*% Sdr.mat    

    ind.mat.minimal.rep <- list()
    for (j in 1:nrow(ind.mat.minimal)) ind.mat.minimal.rep[[j]] <- which.mat(ind.mat.minimal[j,], ind.mat) 
  
    ## break up computation into blocks to not exceed memory limits
    n.per.group <- max(c(round(1e6/d^r),2))
    ##ngroup <- max(n%/%n.per.group+1,1)
    n.seq <- seq(1, n, by=n.per.group)
    if (tail(n.seq,n=1) < n) n.seq <- c(n.seq, n+1)

    if (length(n.seq)> 1)
    {
      mvh.minimal <- numeric()
      for (i in 1:(length(n.seq)-1))
      {
        Hr <- 0
        for (j in 0:floor(r/2))
        {
          cj <- (-1)^j*factorial(r)/(factorial(j)*2^j*factorial(r-2*j))
          vSigmaj <- Kpow(vSigma,j)
          dens <- matrix(exp(logretval[n.seq[i]:(n.seq[i+1]-1)]), nrow=n.seq[i+1]- n.seq[i])
          xr2j <- mat.Kpow(A=x.centred[n.seq[i]:(n.seq[i+1]-1),], pow=r-2*j)
          Hr <- Hr + cj*mat.Kprod(U=xr2j, V=matrix(rep(vSigmaj, n.seq[i+1]- n.seq[i]), nrow=n.seq[i+1]- n.seq[i],byrow=TRUE))
        }
        mvh.temp <- (-1)^r*(t(ones) %x% dens)*Hr %*% SirSdr
        mvh.minimal <- rbind(mvh.minimal, mvh.temp[,ind.mat.minimal.logical])
      }
    }
    else
    {
      Hr <- 0
      for (j in 0:floor(r/2))
      {
        cj <- (-1)^j*factorial(r)/(factorial(j)*2^j*factorial(r-2*j))
        vSigmaj <- Kpow(vSigma,j)
        dens <- matrix(exp(logretval), nrow=n)
        xr2j <- mat.Kpow(A=x.centred, pow=r-2*j)
        Hr <- Hr + cj*mat.Kprod(U=xr2j,V=matrix(rep(vSigmaj,n),nrow=n,byrow=TRUE))
      }
      mvh <- (-1)^r*(t(ones) %x% dens)*Hr %*% SirSdr
      mvh.minimal <- mvh[,ind.mat.minimal.logical]
    }

    if (is.vector(mvh.minimal)) mvh.minimal <- matrix(mvh.minimal, nrow=1)
  }

  
  if (r>0)
  {  
    if (deriv.vec)
    {
      mvh <- matrix(0, nrow=n, ncol=d^r)
      for (j in 1:length(ind.mat.minimal.rep)) mvh[,ind.mat.minimal.rep[[j]]] <- mvh.minimal[,j]
    }
    if (!deriv.vec)
    {
      mvh <- mvh.minimal
      ind.mat <- ind.mat.minimal
    }
  }
  
  if (add.index) return(list(deriv=mvh, deriv.ind=ind.mat))
  else return(deriv=mvh)
}



dmvnorm.deriv.mixt <- function(x, mus, Sigmas, props, deriv.order, Sdr.mat, deriv.vec=TRUE, add.index=FALSE, only.index=FALSE)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  if (is.vector(x)) d <- length(x)
  else d <- ncol(x)
  
  if (missing(mus)) mus <- rep(0,d)
  if (missing(Sigmas)) Sigmas <- diag(d)

  r <- deriv.order
  sumr <- sum(r)
  if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d, r=sumr)
  
  ind.mat <- dmvnorm.deriv(x=x, mu=mus[1,], Sigma=Sigmas[1:d,], deriv.order=r, Sdr.mat=Sdr.mat, only.index=TRUE)
  if (only.index)
    if (deriv.vec) return (ind.mat)
    else return(unique(ind.mat))
  
  ## derivatives 
  ## single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    if (is.matrix(mus)) mus <- mus[1,]
    dens <- dmvnorm.deriv(x=x, mu=mus, Sigma=Sigmas[1:d,], deriv.order=sumr, Sdr.mat=Sdr.mat)
  }
  ## multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0
    ## sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dmvnorm.deriv(x=x, mu=mus[i,], Sigma=Sigmas[((i-1)*d+1):(i*d),], deriv.order=sumr, Sdr.mat=Sdr.mat)  
  }

  if (!deriv.vec)
  { 
    dens <- dens[,!duplicated(ind.mat)]
    ind.mat <- unique(ind.mat)
  }
 
  if (add.index) return(list(deriv=dens, deriv.ind=ind.mat))
  else return(deriv=dens)
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



dmvnorm.deriv.sum <- function(x, Sigma, deriv.order=0, inc=1, binned=FALSE, bin.par, bgridsize, kfe=FALSE, deriv.vec=TRUE, add.index=FALSE, double.loop=FALSE, Sdr.mat, verbose=FALSE)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)
  ind.mat <- dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=diag(d), deriv.order=r, only.index=TRUE)

  if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d,r=r)
  if (missing(bgridsize)) bgridsize <- default.bgridsize(d)
  
  if (binned)
  {
    ##if (!is.diagonal(Sigma)) stop("Binned estimation defined for diagonal Sigma only")
    d <- ncol(Sigma)
    n <- nrow(x)

    if (is.diagonal(Sigma))
    {
      if (missing(bin.par)) bin.par <- binning(x, H=diag(diag(Sigma)), bgridsize=bgridsize)  
      est <- kdde.binned(x=x, bin.par=bin.par, H=Sigma, deriv.order=r, Sdr.mat=Sdr.mat, verbose=verbose)$estimate 
      if (r>0)
      {
        sumval <- rep(0, length(est))
        for (j in 1:length(est)) sumval[j] <- sum(bin.par$counts * n * est[[j]])
      }
      else
        sumval <- sum(bin.par$counts * n * est)
    }
    ## transformation approach from Jose E. Chacon 06/12/2010
    else
    {
      Sigmainv12 <- matrix.sqrt(chol2inv(chol(Sigma)))
      y <- x %*% Sigmainv12
      if (missing(bin.par)) bin.par <- binning(x=y, H=diag(d), bgridsize=bgridsize)  

      est <- kdde.binned(x=y, bin.par=bin.par, H=diag(d), deriv.order=r, Sdr.mat=Sdr.mat, verbose=verbose)$estimate
      if (r>0)
      {
        sumval <- rep(0, length(est))
        for (j in 1:length(est)) sumval[j] <- sum(bin.par$counts * n * est[[j]]) 
      }
      else
        sumval <- sum(bin.par$counts * n * est)

      sumval <- det(Sigmainv12) * sumval  %*% Kpow(Sigmainv12, pow=r)
    }
  }
  else
  {
    if (verbose) pb <- txtProgressBar() 
    if (double.loop)
    {
      ngroup <- n
      sumval <- 0
      for (i in 1:nrow(x))
      {
        if (verbose) setTxtProgressBar(pb, i/ngroup) 
        sumval <- sumval + apply(dmvnorm.deriv(x, mu=x[i,], Sigma=Sigma, deriv.order=r, deriv.vec=deriv.vec, Sdr.mat=Sdr.mat), 2 , sum)
      }
    }
    else
    {
      n.per.group <- max(c(round(1e6/(n*d^r)),1))
      ngroup <- max(n%/%n.per.group+1,1)
      sumval <- 0
      n.seq <- seq(1, n, by=n.per.group)
      if (tail(n.seq,n=1) < n) n.seq <- c(n.seq, n+1)

      if (length(n.seq)> 1)
      {
        for (i in 1:(length(n.seq)-1))
        {  
          difs <- differences(x=x, y=x[n.seq[i]:(n.seq[i+1]-1),])
          sumval <- sumval + apply(dmvnorm.deriv(x=difs, mu=rep(0,d), Sigma=Sigma, deriv.order=r, deriv.vec=deriv.vec, Sdr.mat=Sdr.mat), 2 ,sum)
          if (verbose) setTxtProgressBar(pb, i/ngroup) 
        }
      }
     else
     {
       sumval <- apply(dmvnorm.deriv(x=x, mu=x, Sigma=Sigma, deriv.order=r, deriv.vec=deriv.vec, Sdr.mat=Sdr.mat), 2 , sum)
     }
    }
    if (verbose) close(pb)
  }
  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv(x=rep(0,d), mu=rep(0,d), Sigma=Sigma, deriv.order=r, deriv.vec=deriv.vec, Sdr.mat=Sdr.mat)
  
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
    else sumval <- sumval/(n*(n-1))

  if (add.index)
  {
    if (deriv.vec) return(list(sum=sumval, deriv.ind=ind.mat))
    else return(list(sum=sumval, deriv.ind=unique(ind.mat)))
  }
  else return(sum=sumval)
}

## Single partial derivative of the multivariate normal with scalar variance matrix sigma^2 I_d  
## Code by Jose  Chacon 04/09/2007


dmvnorm.deriv.scalar <- function(x, mu, sigma, deriv.order, binned=FALSE)
{
  r <- deriv.order
  d <- ncol(x)

  sderiv <- sum(r)
  arg <- x/sigma
  darg <- dmvnorm(arg, mean=mu)/(sigma^(sderiv+d))
  for (j in 1:d)
  {
    hmold0 <- 1
    hmold1 <- arg[,j]
    hmnew <- 1
    if (r[j] ==1){hmnew<-hmold1}
    if (r[j] >= 2) ## Multiply by the corresponding Hermite polynomial, coordinate-wise, using Fact C.1.4 in W&J (1995) and Willink (2005, p.273)
      for (i in (2:r[j]))
      {
        hmnew <- arg[,j] * hmold1 - (i - 1) * hmold0
        hmold0 <- hmold1
        hmold1 <- hmnew
      }
    darg <- hmnew * darg
  }
  
  val <- darg*(-1)^sderiv
  return(val)
}


dmvnorm.deriv.scalar.sum <- function(x, sigma, deriv.order=0, inc=1, kfe=FALSE, double.loop=FALSE, binned=FALSE, bin.par)
{
  r <- deriv.order
  d <- ncol(x)
  n <- nrow(x)

  if (binned)
  {
    if (missing(bin.par)) bin.par <- binning(x, H=diag(d)*sigma^2)  
    n <- sum(bin.par$counts)

    ind.mat <- dmvnorm.deriv(x=rep(0,d), Sigma=diag(d), deriv.order=sum(r), deriv.vec=TRUE, only.index=TRUE)
    fhatr <- kdde.binned(bin.par=bin.par, H=sigma^2*diag(d), deriv.order=sum(r), deriv.vec=TRUE, w=rep(1,n), deriv.index=which.mat(r=r, ind.mat)[1])
    ##fhatr <- drvkde(x=bin.par$counts, drv=r, bandwidth=sigma, binned=TRUE, se=FALSE)
    sumval <- sum(bin.par$counts * n * fhatr$est[[1]])
  }
  else
  {
    if (double.loop)
    {
      sumval <- 0
      for (i in 1:nrow(x))
        sumval <- sumval + sum(dmvnorm.deriv.scalar(x, mu=x[i,], sigma=sigma, deriv.order=r))
    }
    else
    {
      ngroup <- round(n^2/1e6)+1
      nn <- n %/% ngroup
      sumval <- 0
      for (i in 1:ngroup)
      {
        difs <- differences(x=x, y=x[((i-1)*nn+1):(i*nn),])
        sumval <- sumval + sum(dmvnorm.deriv.scalar(x=difs, mu=rep(0,d), sigma=sigma, deriv.order=r))
      }
      if (n %% ngroup >0)
      {
        difs <- differences(x=x, y=x[(ngroup*nn+1):n,])
      sumval <- sumval + sum(dmvnorm.deriv.scalar(x=difs, mu=rep(0,d), sigma=sigma, deriv.order=r))
      }
    }
  }

  if (inc==0)
    sumval <- sumval - n*dmvnorm.deriv.scalar(x=t(as.matrix(rep(0,d))), mu=rep(0,d), sigma=sigma, deriv.order=r)
  
  if (kfe)
    if (inc==1) sumval <- sumval/n^2
      else sumval <- sumval/(n*(n-1))
  
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


plotmixt <- function(mus, Sigmas, props, dfs, dist="normal", draw=TRUE, ...)
{
  if (ncol(Sigmas)==2)
    plotmixt.2d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, draw=draw, ...)
  else if (ncol(Sigmas)==3)
    plotmixt.3d(mus=mus, Sigmas=Sigmas, props=props, dfs=dfs, dist=dist, draw=draw, ...) 
}


plotmixt.2d <- function(mus, Sigmas, props, dfs, dist="normal",
    xlim, ylim, gridsize, display="slice", cont=c(25,50,75), abs.cont,
    lty, xlab="x", ylab="y", zlab="Density function",
    theta=-30, phi=40, d=4, add=FALSE, drawlabels=TRUE, nrand=1e5, draw=TRUE, ...)
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
    dens <- dmvnorm.mixt(xy, mus=mus, Sigmas=Sigmas, props=props)
  else if (dist=="t")
    dens <- dmvt.mixt(xy, mus=mus, Sigmas=Sigmas, props=props, dfs=dfs)

  dens.mat <- matrix(dens, ncol=length(x), byrow=FALSE)

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
    hts <- quantile(dens.rand, prob=(100 - cont)/100)
  else
    hts <- abs.cont

  if (draw)
  {  
    disp <- substr(display,1,1)
    
    if (disp=="p")
      persp(x, y, dens.mat, theta=theta, phi=phi, d=d, xlab=xlab, ylab=ylab, zlab=zlab, ...)

    else if (disp=="s")
    {
      if (!add)
        plot(x, y, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
      if (missing(lty))
        lty <- 1
  
      for (i in 1:length(hts)) 
      {
        scale <- cont[i]/hts[i]
        if (missing(abs.cont))
          contour(x, y, dens.mat*scale, level=hts[i]*scale, add=TRUE, drawlabels=drawlabels, lty=lty, ...)
        else
          contour(x, y, dens.mat, level=hts[i], add=TRUE, drawlabels=drawlabels, lty=lty, ...)
      }
    }
    
    else if (disp=="i")
      image(x, y, dens.mat, xlab=xlab, ylab=ylab, ...)
    else if (disp=="f")
    filled.contour(x, y, dens.mat, xlab=xlab, ylab=ylab, ...)
  }
  
  if (exists("hts"))
    fhat <- list(eval.points=list(x, y), estimate=dens.mat, cont=hts)
  else
    fhat <- list(eval.points=list(x, y), estimate=dens.mat)
  
  invisible(fhat)
}

plotmixt.3d <- function(mus, Sigmas, props, dfs, cont=c(25,50,75), abs.cont,
    dist="normal", xlim, ylim, zlim, gridsize, alphavec, colors, add=FALSE, nrand=1e5, draw=TRUE, ...)
{
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
      dens <- dmvnorm.mixt(cbind(xy, z[i]), mus=mus, Sigmas=Sigmas, props=props)
    else if (dist=="t")
      dens <- dmvt.mixt(cbind(xy, z[i]), mus=mus, Sigmas=Sigmas, dfs=dfs, props=props)
    
    dens.mat <- matrix(dens, ncol=length(x), byrow=FALSE)
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
  
  if (draw)
  {
    require(rgl)
    require(misc3d)
 
    
    nc <- length(hts)
    if (missing(colors))
      colors <- rev(heat.colors(nc))
  
    if (missing(alphavec))
      alphavec <- seq(0.1,0.5,length=nc)
    
    if (!add) clear3d()
    
    for (i in 1:nc) 
      contour3d(dens.array, level=hts[nc-i+1], x,y,z, add=TRUE, color=colors[i], alpha=alphavec[i],...)
    decorate3d(...)
    }

  if (exists("hts"))
    fhat <- list(eval.points=list(x, y, z), estimate=dens.array, cont=hts)
  else
    fhat <- list(eval.points=list(x, y, z), estimate=dens.array)
  
  invisible(fhat)
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
    dens <- dmvt(x, delta=mus, sigma=Sigmas, df=dfs, log=FALSE)
  
  ## multiple component mixture
  else   
  {   
    if (is.vector(mus)) d <- length(mus)
    else d <- ncol(mus)
    k <- length(props)
    dens <- 0      
    for (i in 1:k)
      dens <- dens+props[i]*dmvt(x,delta=mus[i,],sigma=Sigmas[((i-1)*d+1):(i*d),], df=dfs[i], log=FALSE)
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


