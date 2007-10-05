
###############################################################################
# Univariate mixture normal densities
###############################################################################


rnorm.mixt <- function(n=100, mus=0, sigmas=1, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    rand <- rnorm(n=n, mean=mus, sd=sigmas)

  ### multiple component mixture
  else
  {
    k <- length(props)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    # compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k)
    {
      # compute random sample from normal mixture component
      if (n.prop[i] > 0)
          rand <- c(rand, rnorm(n=n.prop[i], mean=mus[i], sd=sigmas[i]))        
    }
  }

  return(rand[sample(n)])
}


dnorm.mixt <- function(x, mus=0, sigmas=1, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  # single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm(x, mean=mus, sd=sigmas)

  # multiple component mixture
  else   
  {   
    k <- length(props)
    dens <- 0

    # sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dnorm(x, mean=mus[i], sd=sigmas[i])
  }
  
  return(dens)
}   


###############################################################################
# Multivariate normal densities and derivatives
###############################################################################


###############################################################################
# Multivariate normal mixture - random sample
# 
# Parameters
# n - number of samples
# mus - matrix of means (each row is a vector of means from each component
#       density)
# Sigmas - matrix of covariance matrices (every d rows is a covariance matrix 
#          from each component density) 
# props - vector of mixing proportions 
# 
# Returns
# Vector of n observations from the normal mixture 
###############################################################################

rmvnorm.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  #if (is.vector(Sigmas))
  #  return(rnorm.mixt(n=n, mus=mus, sigmas=Sigmas, props=props))
  
  ### single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    rand <- rmvnorm(n=n, mean=mus, sigma=Sigmas)

  ### multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    # compute number taken from each mixture
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))
    
    rand <- numeric(0)
    
    for (i in 1:k)
    {
      # compute random sample from normal mixture component
      if (n.prop[i] > 0)
          rand <- rbind(rand, rmvnorm(n=n.prop[i], mean=mus[i,], 
                                      sigma=Sigmas[((i-1)*d+1) : (i*d),]))        
    }
  }

  return(rand[sample(n),])
}


###############################################################################
# Multivariate normal mixture - density values
# 
# Parameters
# x - points to compute density at 
# mus - matrix of means
# Sigmas - matrix of covariance matrices 
# props - vector of mixing proportions 
# 
# Returns
# Density values from the normal mixture (at x)
###############################################################################

dmvnorm.mixt <- function(x, mus, Sigmas, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")

  #if (is.vector(Sigmas))
  #  return(dnorm.mixt(x=x, mus=mus, sigmas=Sigmas, props=props))
  
  # single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dmvnorm(x, mean=mus, sigma=Sigmas)

  # multiple component mixture
  else   
  {   
    if (is.vector(mus)) d <- length(mus)
    else d <- ncol(mus)
    k <- length(props)
    dens <- 0

    # sum of each normal density value from each component at x  
    for (i in 1:k)
      dens <- dens + props[i]*dmvnorm(x, mean=mus[i,],
                                      sigma=Sigmas[((i-1)*d+1):(i*d),])
  }
  
  return(dens)
}   

###############################################################################
# Compute moments of multivariate normal mixture
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
# Double sum  of K(X_i - X_j) used in density derivative estimation
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.1d.sum <- function(x, sigma, inc=1, binned=FALSE, bin.par)
{
  n <- length(x)
  d <- 1
  if (binned)
  {
    fhatr <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=sigma,
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.1d(x=rep(0,d), r=rep(0,d), Sigma)
  }
  else
  {  
    sumval <- 0
    for (i in 1:n)
      sumval <- sumval + sum(dnorm(x[i] - x, mean=0, sd=sigma))
    
    if (inc == 0) 
      sumval <- sumval - n*dnorm(0, mean=0, sd=sigma)
   } 
  
  return(sumval)
}

###############################################################################
# Double sum  of K(X_i - X_j) used in density derivative estimation
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.2d.sum <- function(x, Sigma, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; d <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); d <- 2; x1 <- x[,1]; x2 <- x[,2]
  }

  if (binned)
  {
    fhatr <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.2d(x=rep(0,d), r=rep(0,d), Sigma)
  }
  else
  { 
    viSigma <- vec(chol2inv(chol(Sigma)))
    result <- .C("dmvnorm_2d_sum", as.double(x1), as.double(x2),
                 as.double(viSigma), as.double(det(Sigma)), as.integer(n),
                 as.double(0), PACKAGE="ks")
    sumval <- result[[6]]
    
    ## above C function mvnorm_2d_sum only computes the upper triangular half
    ## so need to reflect along the diagonal and then subtract appropriate
    ## amount to compute whole sum 
  
    if (inc == 0) 
      sumval <- 2*sumval - 2*n*dmvnorm(c(0,0), c(0,0), Sigma)
    else if (inc == 1)
      sumval <- 2*sumval - n*dmvnorm(c(0,0), c(0,0), Sigma) 
  }
  
  return(sumval)
}

###############################################################################
# Double sum (i = 1:n1, j = 1:n2) of K(X_i - X_j), X_i from cluster 1,
# X_j from cluster 2
#
# Parameters
# x - data points from first cluster
# y - data points from second cluster
# Sigma - variance matrix
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.2d.sum.pc <- function(x, y, Sigma, inc=1)
{
  if (is.vector(x))
  {
    n1 <- 1; d <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n1 <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
  if (is.vector(y))
  {
    n2 <- 1; d <- 1; y1 <- y[1]; y2 <- y[2]
  }
  else
  {
    n2 <- nrow(y); y1 <- y[,1]; y2 <- y[,2]
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_2d_sum_clust", as.double(x1), as.double(x2),
               as.double(y1), as.double(y2), as.double(viSigma),
               as.double(det(Sigma)), as.integer(n1), as.integer(n2),
               as.double(0), PACKAGE="ks")
  sumval <- result[[9]]
  
  return(sumval)
}  

###############################################################################
# Partial derivatives of the univariate normal (mean 0) 
# 
# Parameters
# x - points to evaluate at
# sigma - std deviation
# r - derivative index 
#
# Returns
# r-th derivative at x
###############################################################################

dmvnorm.deriv.1d <- function(x, sigma, r)
{

  phi <- dnorm(x, mean=0, sd=sigma) 
  if (r==0)
    return(phi)
  else if (r==1)
    derivt <- x*phi
  else if (r==2)
    derivt <- (x^2-1)*phi
  else if (r==3)
    derivt <- (x^3 - x)*phi
  else if (r==4)
    derivt <- (x^4 - 6*x^2 + 3)*phi
  else if (r==5)
    derivt <- (x^5 - 10*x^3 + 15*x)*phi
  else if (r==6)
    derivt <- (x^6 - 15*x^4 + 45*x^2 -15)*phi
  else if (r==7)
    derivt <- (x^7 - 21*x^5 + 105*x^3 - 105*x)*phi
  else if (r==8)
    derivt <- (x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105)*phi
  else
    stop("Only works for up to 8th order partial derivatives")
  
  return(derivt)  
}

###############################################################################
# Partial derivatives of the bivariate normal (mean 0) 
# 
# Parameters
# x - points to evaluate at
# Sigma - variance matrix
# r - (r1, r2) vector of partial derivative indices 
#
# Returns
# r-th partial derivative at x
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
  # first order derivatives  
  else if (sum(r) == 1)
    derivt <- .C("dmvnormd1_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  # third order derivatives
  else if (sum(r) == 2)
    derivt <- .C("dmvnormd2_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  # second order derivatives
  else if (sum(r) == 3)
    derivt <- .C("dmvnormd3_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  # fourth order derivatives
  else if (sum(r) == 4)
    derivt <- .C("dmvnormd4_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  # fifth order derivatives
  else if (sum(r) == 5)
    derivt <- .C("dmvnormd5_2d", as.double(x1), as.double(x2), 
                 as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  # sixth order derivatives
  else if (sum(r) == 6)
    derivt <- .C("dmvnormd6_2d", as.double(x1), as.double(x2), 
            as.double(vec(Sigma)), as.integer(r), as.integer(n), 
                 as.double(rep(0, n)), PACKAGE="ks")
  else
    stop("Only works for up to 6th order partial derivatives")
  
  return(derivt[[6]])
}         



###############################################################################
# Double sum of K^(r)(X_i - X_j) used in density derivative estimation
#
# Parameters
# x - points to evaluate 
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.deriv.2d.sum <- function(x, Sigma, r, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
  
  if (binned)
  {
    ## drvkde computes include-diagonals estimate
    fhatr <- drvkde.ks(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma)
  }
  else  
  {
    if (sum(r)==0)
    return(dmvnorm.2d.sum(x, Sigma, inc=inc))
    ## 1st order derivatives
    else if (sum(r)==1)
      derivt <- .C("dmvnormd1_2d_sum", as.double(x1), as.double(x2), 
                   as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0),
                   PACKAGE="ks")
    ## 2nd order derivativeselse if (sum(r)==2)
    else if (sum(r)==2)
      derivt <- .C("dmvnormd2_2d_sum", as.double(x1), as.double(x2), 
                   as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0),
                   PACKAGE="ks")
    ## 3rd order derivativeselse if (sum(r)==3)
    else if (sum(r)==3)
      derivt <- .C("dmvnormd3_2d_sum", as.double(x1), as.double(x2), 
                   as.double(vec(Sigma)), as.integer(r), as.integer(n), as.double(0),
                   PACKAGE="ks")
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
      sumval <- 2*sumval - 2*n*dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma)
    else if (inc == 1)
      sumval <- 2*sumval - n*dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma) 
  }
      
  return(sumval)
}         



###############################################################################
# Double sum K^(r)(X_i - X_j)*(X_i - X_j)*(X_i - X_j)^T used in
# density derivative estimation (Hbcv)
#
# Parameters
# x - points to evaluate 
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.deriv.2d.xxt.sum <- function(x, Sigma, r)
{
  if (is.vector(x))
  {
    n <- 1; d <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnormd4_2d_xxt_sum", as.double(x1), as.double(x2),
               as.double(viSigma), as.integer(r), as.integer(n), 
               as.double(rep(0,4)), PACKAGE="ks")
  
  # above C functions dmvnorm4_2d_xxt_sum only computes the upper triangular
  # half so need to reflect along the diagonal to compute whole sum
  sumval <- 2*result[[6]]
  
  return(invvec(sumval))
} 
###############################################################################
# Double sum  of K(X_i - X_j) used in density derivative estimation - 3-dim
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.3d.sum <- function(x, Sigma, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3];
  }

  if (binned)
  {
     fhatr <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=sqrt(diag(Sigma)),
                        binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
     sumval <- sum(bin.par$counts * n * fhatr)
     if (inc == 0) 
       sumval <- sumval - n*dmvnorm.deriv.3d(x=rep(0,d), r=rep(0,d), Sigma)
  }
  else
  {  
    viSigma <- vec(chol2inv(chol(Sigma)))
    result <- .C("dmvnorm_3d_sum", as.double(x1), as.double(x2), as.double(x3), 
                 as.double(viSigma), as.double(det(Sigma)), as.integer(n),
                 as.double(0), PACKAGE="ks")
    sumval <- result[[7]]

    ## above C function mvnorm_3d_sum only computes the upper triangular half
    ## so need to reflect along the diagonal and then subtract appropriate
    ## amount to compute whole sum 
  
    if (inc == 0) 
      sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
    else if (inc == 1)
      sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  }
  
  return(sumval)
}




###############################################################################
# Partial derivatives of the 3-variate normal (mean 0, diag variance matrix) 
# 
# Parameters
# x - points to evaluate at
# Sigma - variance matrix
# r - (r1, r2, r3, r4) vector of partial derivative indices 
#
# Returns
# r-th partial derivative at x
##############################################################################

dmvnorm.deriv.3d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment
  
  d <- 3
  if (sum(r) > 6)
    stop("Only works for up to 6th order partial derivatives")

  if (is.vector(x))
  {    
    n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; 
  }
  else 
  {
    n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; 
  }
  
  y1 <- cbind(x1,x2)
  y2 <- x3
  r1 <- r[c(1,2)]
  r2 <- r[3]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- Sigma[3,3]

  # Use existing 2-dim derivatives to compute 3-dim derivative for diag
  # variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.1d(y2, sqrt(Sigma2), r2)
  derivt <- derivt1*derivt2

  return(derivt)
}         

###############################################################################
# Double sum of K^(r)(X_i - X_j) used in density derivative estimation - 3-dim
#
# Parameters
# x - points to evaluate 
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# Double sum at x
###############################################################################

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
    fhatr <- drvkde.ks(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.3d(x=rep(0,d), r=r, Sigma)
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
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

dmvnorm.4d.sum <- function(x, Sigma, inc=1, binned=FALSE, bin.par)
{
  if (is.vector(x))
  {
    n <- 1; d <- 4; x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; 
  }
  else
  {
    n <- nrow(x); d <- ncol(x); x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
  }

  if (binned)
  {
     fhatr <- drvkde.ks(x=bin.par$counts, drv=rep(0,d), bandwidth=sqrt(diag(Sigma)),
                        binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
     sumval <- sum(bin.par$counts * n * fhatr)
     if (inc == 0) 
       sumval <- sumval - n*dmvnorm.deriv.4d(x=rep(0,d), r=rep(0,d), Sigma)
  }
  else
  {
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
      sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
    else if (inc == 1)
      sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  }

  return(sumval)
}


###############################################################################
# Partial derivatives of the 4-variate normal (mean 0, diag variance matrix) 
# 
# Parameters
# x - points to evaluate at
# Sigma - variance matrix
# r - (r1, r2, r3, r4) vector of partial derivative indices 
#
# Returns
# r-th partial derivative at x
##############################################################################

dmvnorm.deriv.4d <- function(x, Sigma, r) 
{
  ####### Diagonal variance matrix implemented ONLY at the moment
  
  d <- 4
  if (sum(r) > 6)
    stop("Only works for up to 6th order partial derivatives")

  if (is.vector(x))
  {    
    n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
  }
  else 
  {
    n <- nrow(x);
    x1 <- x[,1]; x2 <- x[,2]; x3 <- x[,3]; x4 <- x[,4];
  }
  
  y1 <- cbind(x1,x2)
  y2 <- cbind(x3,x4)
  r1 <- r[c(1,2)]
  r2 <- r[c(3,4)]
  Sigma1 <- diag(c(Sigma[1,1], Sigma[2,2]))
  Sigma2 <- diag(c(Sigma[3,3], Sigma[4,4]))

  # Use existing 2-dim derivatives to compute 4-dim derivative for diag
  # variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt <- derivt1*derivt2

  return(derivt)
}         

###############################################################################
# Double sum of K^(r)(X_i - X_j) used in density derivative estimation - 4-dim
#
# Parameters
# x - points to evaluate 
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# Double sum at x
###############################################################################

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
    fhatr <- drvkde.ks(x=bin.par$counts, drv=r, bandwidth=sqrt(diag(Sigma)),
                       binned=TRUE, range.x=bin.par$range.x, se=FALSE)$est
    sumval <- sum(bin.par$counts * n * fhatr)
    if (inc == 0) 
      sumval <- sumval - n*dmvnorm.deriv.4d(x=rep(0,d), r=r, Sigma)
  }
  else
  {  
    for (j in 1:n)
    {  
      y1 <- x1 - x1[j]
      y2 <- x2 - x2[j]
      y3 <- x3 - x3[j]
      y4 <- x4 - x4[j]
      sumval <- sumval + sum(dmvnorm.deriv.4d(cbind(y1, y2, y3, y4), Sigma, r))
    }
  
    if (inc==0)
      sumval <- sumval - n*dmvnorm.deriv.4d(rep(0,d), Sigma, r)
  }
  
  return(sumval)
}


###############################################################################
# Double sum  of K(X_i - X_j) used in density derivative estimation - 5-dim
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
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

  # above C function mvnorm_5d_sum only computes the upper triangular half
  # so need to reflect along the diagonal and then subtract appropriate
  # amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  
  return(sumval)
}

###############################################################################
# Partial derivatives of the 5-variate normal (mean 0, diag variance matrix) 
# 
# Parameters
# x - points to evaluate at
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# r-th partial derivative at x
##############################################################################

dmvnorm.deriv.5d <- function(x, Sigma, r) 
{
  d <- 5
  if (sum(r) > 6)
    stop("Only works for 2nd, 4th and 6th order partial derivatives")
 
  if (is.vector(x))
  {    
    n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]; 
  }
  else 
  {
    n <- nrow(x);
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
  
  # Use existing 2-dim derivatives to compute 6-dim derivative for diag
  # variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt3 <- dmvnorm.deriv.1d(y3, sqrt(Sigma3), r3)
  derivt <- derivt1*derivt2*derivt3

  return(derivt)
}         


###############################################################################
# Double sum  of K^(r)(X_i - X_j) used in density derivative estimation - 5-dim
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

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
# Double sum  of K(X_i - X_j) used in density derivative estimation - 6-dim
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
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
    x5 <- x[,5]; x6 <- x[6];  
  }
  
  viSigma <- vec(chol2inv(chol(Sigma)))
  result <- .C("dmvnorm_6d_sum", as.double(x1), as.double(x2),
               as.double(x3), as.double(x4), as.double(x5), as.double(x6),
               as.double(viSigma), as.double(det(Sigma)), as.integer(n),
               as.double(0), PACKAGE="ks")
  sumval <- result[[10]]

  # above C function mvnorm_6d_sum only computes the upper triangular half
  # so need to reflect along the diagonal and then subtract appropriate
  # amount to compute whole sum 
  
  if (inc == 0) 
    sumval <- 2*sumval - 2*n*dmvnorm(rep(0,d), rep(0,d), Sigma)
  else if (inc == 1)
    sumval <- 2*sumval - n*dmvnorm(rep(0,d), rep(0,d), Sigma) 
  
  return(sumval)
}


###############################################################################
# Partial derivatives of the 6-variate normal (mean 0, diag variance matrix) 
# 
# Parameters
# x - points to evaluate at
# Sigma - variance matrix
# r - vector of partial derivative indices 
#
# Returns
# r-th partial derivative at x
##############################################################################

dmvnorm.deriv.6d <- function(x, Sigma, r) 
{
  d <- 6
  if (sum(r) > 6)
    stop("Only works for 2nd, 4th and 6th order partial derivatives")

  if (is.vector(x))
  {    
    n <- 1;
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]; x5 <- x[5]; x6 <- x[6];
  }
  else 
  {
    n <- nrow(x);
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
  
  # Use existing 2-dim derivatives to compute 6-dim derivative for diag
  # variance matrix
  derivt1 <- dmvnorm.deriv.2d(y1, Sigma1, r1)
  derivt2 <- dmvnorm.deriv.2d(y2, Sigma2, r2)
  derivt3 <- dmvnorm.deriv.2d(y3, Sigma3, r3)
  derivt <- derivt1*derivt2*derivt3

  return(derivt)
}         


###############################################################################
# Double sum  of K^(r)(X_i - X_j) used in density derivative estimation - 6-dim
#
# Parameters
# x - points to evaluate
# Sigma - variance matrix
# inc - 0 - exclude diagonals
#     - 1 - include diagonals
#
# Returns
# Double sum at x
###############################################################################

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
# Multivariate t - density values
#
# Parameters
# x - points to compute density     
# mu - vector of means 
# Sigma - dispersion matrix
# df - degrees of freedom
#
# Returns
# Value of multivariate t density at x
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
# Multivariate t mixture - density values
#
# Parameters
# x - points to compute density at    
# mus - vector of means 
# Sigmas - dispersion matrices
# dfs - degrees of freedom
# props - vector of mixing proportions
#
# Returns
# Value of multivariate t mixture density at x
###############################################################################

dmvt.mixt <- function(x, mus, Sigmas, dfs, props)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")
  
  # single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dmvt(x, mu=mus, Sigma=Sigmas, df=dfs)
  
  # multiple component mixture
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
# Multivariate t mixture - random sample
# 
# Parameters
# n - number of samples
# mus - means 
# Sigmas - matrix of dispersion matrices
# dfs - vector of degrees of freedom
# props - vector of mixing proportions 
# 
# Returns
# Vector of n observations from the t mixture
###############################################################################

rmvt.mixt <- function(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))  
    stop("Proportions don't sum to one\n")
  else if (length(dfs) != length(props))
    stop("Length of df and mixing proportions vectors not equal")  

  # single component mixture
  if (identical(all.equal(props[1], 1), TRUE))
  {
    rand <- rmvt(n=n, sigma=Sigmas, df=dfs)
    for (i in 1:length(mus))
      rand[,i] <- rand[,i] + mus[i]
  }
  
  # multiple component mixture
  else
  {
    k <- length(props)
    d <- ncol(Sigmas)
    n.samp <- sample(1:k, n, replace=TRUE, prob=props) 
    n.prop <- numeric(0)

    # compute number to be drawn from each component 
    for (i in 1:k)
      n.prop <- c(n.prop, sum(n.samp == i))

    # generate random samples from each component
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


###############################################################################
# Creates plots of mixture density functions
#
# Parameters
# mus - means
# Sigmas - variances
# props - vector of proportions of each mixture component 
# dfs - degrees of freedom
# dist - "normal" - normal mixture
#      - "t" - t mixture
# ...
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
    scale <- cont[i]/hts[i]
    contour3d(dens.array, level=hts[nc-i+1],x, y, z, add=TRUE, color=colors[i],
             alpha=alphavec[i], ...)
  }
}



