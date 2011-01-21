
#############################################################################
## Basic important vectors and matrices
#############################################################################
  
###############################################################################
# Vec operator 
#
# Takes square matrix x and returns its elements stacked column wise in a
# vector
###############################################################################

vec <- function(x, byrow=FALSE)
{
  if (is.vector(x)) return (x)
  
  if (byrow) x <- t(x)
  d <- ncol(x)
  vecx <- vector()
  for (j in 1:d)
    vecx <- c(vecx, x[,j])
  
  return(vecx)           
}

###############################################################################
# Vech operator 
# 
# Takes matrix x and returns its elements that in the lower triangular half, 
# stacked column wise in a vector
###############################################################################

vech <- function(x)
{
  if (is.vector(x))
  {
    if (length(x)==1)
      return (x)
    else
      stop("vech undefined for vectors")
  }
  else if (is.matrix(x))
  {  
    d <- ncol(x)
    if (d!=nrow(x)) ##if (!isSymmetric(x))
      stop("vech only defined for square matrices")
    
    vechx <- vector()
    for (j in 1:d)
      vechx <- c(vechx, x[j:d,j])
    return(vechx)
  }
  
}

  
###############################################################################
# Inverse vec operator 
#
# Takes vector x and stacks its elements into columns of a matrix. 
###############################################################################
    
invvec <- function(x, ncol, nrow, byrow=FALSE)
{
  if (length(x)==1)
    return(x)
  
  d <- sqrt(length(x))
  if (missing(ncol) | missing(nrow))
  {
    ncol <- d; nrow <- d
    if (round(d) != d)
      stop("Need to specify nrow and ncol for non-square matrices")
  }
  
  invvecx <- matrix(0, nrow = nrow, ncol = ncol)
  if (byrow)
    for (j in 1:nrow)
      invvecx[j,] <- x[c(1:ncol) + (j-1)*ncol]
  else
    for (j in 1:ncol)
      invvecx[,j] <- x[c(1:nrow) + (j-1)*nrow]
  
  return(invvecx)
}

###############################################################################
# Inverse vech operator 
#
# Takes vector x and stacks its elements into a symmetric matrix 
###############################################################################

invvech <- function(x)
{
  if (length(x)==1)
    return(x)
  
  d <- (-1 + sqrt(8*length(x) + 1))/2
  if (round(d) != d)
    stop("Number of elements in x will not form a square matrix.")
  invvechx <- matrix(0, nrow=d, ncol=d)

  for (j in 1:d)
    invvechx[j:d,j] <- x[1:(d-j+1)+ (j-1)*(d - 1/2*(j-2))]
  
  invvechx <- invvechx + t(invvechx) - diag(diag(invvechx))
  
  return(invvechx)
}

##### Trace of matrix

tr <- function(A)
{
  count <- 0
  if (is.vector(A)) return (A[1])
  if (nrow(A)!=ncol(A))
    stop('Not square matrix')

  else 
     for (i in 1:nrow(A))
       count <- count + A[i,i]
  return(count)
}


###############################################################################
# Elementary vector 
# 
# Creates a vector of length d and with 1 at i-th component and 0 elsewhere 
###############################################################################
    
elem <- function(i, d)
{
  elem.vec <- rep(0, d)
  elem.vec[i] <- 1
  
  return(elem.vec)
}      

### Commutation matrix (taken from MCMCglmmm library)

comm <- function(m,n){
  K<-matrix(0,m*n, m*n)
  H<-matrix(0,m,n)
  for(i in 1:m){
    for(j in 1:n){ 
      H[i,j]<-1
      K<-K+kronecker(H,t(H))
      H[i,j]<-0
    }
  }
  return(K)
}

###############################################################################
# Duplication matrix
# Taken from Felipe Osorio http://www.ime.usp.br/~osorio/files/dupl.q
###############################################################################

dupl <- function(order, ret.q = FALSE)
{
    # call
    cl <- match.call()
    time1 <- proc.time()
    if (!is.integer(order))
        order <- as.integer(order)
    n <- order - 1
    
    # initial duplication matrix
    d1 <- matrix(0, nrow = 1, ncol = 1)
    d1[1,1] <- 1
    if (!is.integer(d1))
        storage.mode(d1) <- "integer"
    
    # recursive formula
    if (n > 0){
    	for (k in 1:n){
    	    drow <- 2*k + 1 + nrow(d1)
    	    dcol <- k + 1 + ncol(d1)
    	    d2 <- matrix(0, nrow = drow, ncol=dcol)
    	    storage.mode(d2) <- "integer"
    	    d2[1,1] <- 1
    	    d2[2:(k+1),2:(k+1)] <- diag(k)
    	    d2[(k+2):(2*k+1),2:(k+1)] <- diag(k)
    	    d2[(2*k+2):drow,(k+2):dcol] <- d1
    	    # permutation matrix
    	    q <- permute.mat(k)
    	    # new duplication matrix
    	    d2 <- q %*% d2
    	    storage.mode(d2) <- "integer"
    	    d1 <- d2
    	}
    }
    else {
    	d2 <- q <- d1
    }
    
    # results
    obj <- list(call=cl, order=order, d=d2)
    if (ret.q)
        obj$q <- q
    obj$time <- proc.time() - time1
    obj
}

invdupl <- function(order, ret.q = FALSE)
{
    # call
    cl <- match.call()
    time1 <- proc.time()
    if (!is.integer(order))
        order <- as.integer(order)
    n <- order - 1

    # initial inverse of duplication matrix
    h1 <- matrix(0, nrow = 1, ncol = 1)
    h1[1,1] <- 1

    # recursive formula
    if (n > 0){
    	for (k in 1:n){
    	    hrow <- k + 1 + nrow(h1)
    	    hcol <- 2*k + 1 + ncol(h1)
    	    h2 <- matrix(0, nrow = hrow, ncol=hcol)
    	    h2[1,1] <- 1
    	    h2[2:(k+1),2:(k+1)] <- .5*diag(k)
    	    h2[2:(k+1),(k+2):(2*k+1)] <- .5*diag(k)
    	    h2[(k+2):hrow,(2*k+2):hcol] <- h1
    	    # permutation matrix
    	    q <- permute.mat(k)
    	    # new inverse of duplication matrix
    	    h2 <- h2 %*% t(q)
    	    h1 <- h2
    	}
    }
    else {
    	h2 <- q <- h1
    }
    
    # results
    obj <- list(call=cl, order=order, h=h2)
    if (ret.q)
        obj$q <- q
    obj$time <- proc.time() - time1
    obj
}




###############################################################################
# Matrix square root - taken from Stephen Lake 
# http://www5.biostat.wustl.edu/s-news/s-news-archive/200109/msg00067.html
###############################################################################

matrix.sqrt <- function(A)
{
  if (length(A)==1)
    return(sqrt(A))
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

matrix.pow <- function(A, n)
{
  if (nrow(A)!=ncol(A)) stop("A must a a square matrix")
  if (floor(n)!=n) stop("n must be an integer")
  if (n==0) return(diag(ncol(A)))
  if (n < 0) return(matrix.pow(A=chol2inv(chol(A)), n=-n))
        
  # trap non-integer n and return an error
  if (n == 1) return(A)
  result <- diag(1, ncol(A))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% A
      n <- n - 1
    }
    A <- A %*% A
    n <- n / 2
  }
  return(result)
}


###############################################################################
# Pre-sphering
# Parameters
# x - data points
#
# Returns
# Pre-sphered x values
###############################################################################

pre.sphere <- function(x, mean.centred=FALSE)
{
  S <- var(x)
  Sinv12 <- matrix.sqrt(chol2inv(chol(S)))

  if (mean.centred)
  {
    xmean <- apply(x,2,mean)
    for (i in 1:ncol(x))
      x[,i] <- x[,i] - xmean[i]
  }
  x.sphered <- x %*% Sinv12
  ##x.sphered <- matrix(0, ncol=ncol(x), nrow=nrow(x))
  ##for (i in 1:nrow(x))
  ##  x.sphered[i,] <- Sinv12 %*% x[i,]    

  return (x.sphered)
}


###############################################################################
# Pre-scaling
# Parameters
# x - data points
#
# Returns
# Pre-scaled x values
###############################################################################

pre.scale <- function(x, mean.centred=FALSE)
{
  x.scaled <- numeric()
  x.sd <- apply(x, 2, sd)
  d <- ncol(x)

  for (i in 1:d)
    if (mean.centred)
      x.scaled <- cbind(x.scaled, (x[,i] - mean(x[,i]))/x.sd[i])
    else
      x.scaled <- cbind(x.scaled, x[,i]/x.sd[i])                  

  return (x.scaled)
}


###############################################################################
# Finds row index matrix
# Parameters
# x - data points
#
# Returns
# i  - if r==mat[i,]
# NA - otherwise
###############################################################################

which.mat <- function(r, mat)
{
  ind <- numeric()
  
  for (i in 1:nrow(mat))
    if (identical(r, mat[i,])) ind <- c(ind,i)

  return(ind)  
}


###############################################################################
# Permute a list of values
#
# Same function as EXPAND.GRID (base package), modified to take 
# list as an argument and returns a matrix 
###############################################################################

permute <- function (args) 
{
  nargs <- length(args)
  if (!nargs) 
    return(as.data.frame(list()))
  if (nargs == 1 && is.list(a1 <- args[[1]])) 
    nargs <- length(args <- a1)
  if (nargs <= 1) 
    return(as.data.frame(if (nargs == 0 || is.null(args[[1]])) list() else args, 
                         optional = TRUE))
  cargs <- args
  rep.fac <- 1
  ##orep <- final.len <- prod(sapply(args, length))
  orep <- prod(sapply(args, length))
  
  for (i in 1:nargs) {
    x <- args[[i]]
    nx <- length(x)
    orep <- orep/nx
    cargs[[i]] <- rep(rep(x, rep(rep.fac, nx)), orep)
    rep.fac <- rep.fac * nx
  }
  do.call("cbind", cargs)
} 

permute.mat <- function(order)
{
    m <- as.integer(order)
    m <- m + 1
    eye <- diag(m)
    u1 <- eye[1:m,1]
    u2 <- eye[1:m,2:m]
    q1 <- kronecker(eye, u1)
    q2 <- kronecker(eye, u2)
    q <- matrix(c(q1, q2), nrow = nrow(q2), ncol = ncol(q1) + ncol(q2))
    if (!is.integer(q))
        storage.mode(q) <- "integer"
    q
}



###############################################################################
## Boolean functions
###############################################################################

is.even <- function(x)
{
  y <- x[x>0] %%2
  return(identical(y, rep(0, length(y))))
}

is.diagonal <- function(x)
{
  return(identical(diag(diag(x)),x))
}


####################################################################
## Functions for unconstrained pilot selectors, and (A)MISE-optimal
## selectors for normal mixtures
## Author: Jose E. Chacon
####################################################################


differences <- function(x, y, upper=FALSE)
{
  if (missing(y)) y <- x
  if (is.vector(x)) x <- t(as.matrix(x))
  if (is.vector(y)) y <- t(as.matrix(y))

  nx <- nrow(x)
  ny <- nrow(y)
  d <- ncol(x)

  difs <- matrix(ncol=d,nrow=nx*ny)
  for (j in 1:d)
  {    
    difs[,j] <- rep(x[,j], times=ny) - rep(y[1:ny,j], each=nx)
    ##difs[,j] <- as.vector(x[,j]%*%t(rep(1,ny))-rep(1,nx)%*%t(y[,j]))
    ##The jth column of difs contains all the differences X_{ij}-Y_{kj}
  }
 
  if (upper)
  {
    ind.remove <- numeric()
    for (j in 1:(nx-1))
      ind.remove <- c(ind.remove, (j*nx+1):(j*nx+j))
      
    return(difs[-ind.remove,])
  }
  else
    return(difs)
}


##### Odd factorial

OF<-function(m){factorial(m)/(2^(m/2)*factorial(m/2))} 


##### Commutation matrix of order m,n

K.mat<-function(m,n){
  K<-0
  for(i in 1:m){for(j in 1:n){
    H<-matrix(0,nrow=m,ncol=n)
    H[i,j]<-1
    K<-K+(H%x%t(H))
  }}
  return(K)        
}

## Kronecker sum

Ksum <- function(A,B)
{
  AB <- numeric()
  for (i in 1:nrow(A))
    for (j in 1:nrow(B))
      AB <- rbind(AB, A[i,] + B[j,])

  return(AB)
}



##### Kronecker power of a matrix A

Kpow<-function(A,pow){    #### WARNING! Dot omitted from original K.pow name
  if(floor(pow)!=pow)stop("pow must be an integer")
  Apow<-A
  if(pow==0){Apow<-1}
  if(pow>1){
    for(i in 2:pow) Apow<-Apow%x%A
  }
  return(Apow)
}
    
##### Row-wise Kronecker products and powers of matrices
    
mat.Kprod<-function(U,V){ #### Returns a matrix with rows U[i,]%x%V[i,]
 n1<-nrow(U)
 n2<-nrow(V)
 if(n1!=n2)stop("U and V must have the same number of vectors")
 p<-ncol(U)
 q<-ncol(V)
 onep<-rep(1,p)
 oneq<-rep(1,q)
 P<-(U%x%t(oneq))*(t(onep)%x%V)
 return(P)
} 


#### WARNING! Dot omitted from original mat.K.pow name
mat.Kpow<-function(A,pow){ #### Returns a matrix with the pow-th Kronecker power of A[i,] in the i-th row
  Apow<-A
  if(pow==0){Apow<-matrix(1,nrow=nrow(A), ncol=1)}  
  if(pow>1){
    for(i in 2:pow) Apow<-mat.Kprod(Apow,A)
  }
  return(Apow)
}


##### Vector of all r-th partial derivatives of the normal density at x=0, i.e., D^{\otimes r)\phi(0)

DrL0<-function(d,r,Sdr){
  DL0<-(-1)^(r/2)*(2*pi)^(-d/2)*OF(r)*(Sdr%*%Kpow(A=vec(diag(d)),pow=r/2))
  return(DL0)
}



T<-function(d,r){    #### Second version, recursive
  Id<-diag(d)
  Tmat<-Id
  Kdd<-K.mat(d,d)
  if(r>1){for(j in 2:r){
    Idj2Kdd<-diag(d^(j-2))%x%Kdd
    Tmat<-Idj2Kdd%*%((Tmat%x%Id)+Idj2Kdd)%*%Idj2Kdd
  }}
  return(Tmat)
}


Trow <- function(d,r,row=1){    #### Second version, recursive
  Id <- diag(d)
  Tmat <- Id
  Kdd <- K.mat(d,d)
  if (r>1)
  {
    for (j in 2:r)
    {
     Idj2Kdd <- diag(d^(j-2)) %x% Kdd
     browser()
     Tmat <- Idj2Kdd %*% ((Tmat %x% Id) + Idj2Kdd) %*% Idj2Kdd
   }
  } 
  return(Tmat)
}

## symmetriser matrix

Sdr<-function(d, r){
  if (r==0) S<-1
  else{
    S <- diag(d)
    if (r>=2)
      for(j in 2:r)
        S <- S %x% diag(d) %*% T(d,j)/j
  }
  return(S)
}

Sdrx <- function(d,r,x)
{
  xtemp <- x
    
}

########################################################################
### Identifying elements of Theta_6 matrix
########################################################################

Theta6.elem <- function(d)
{
  Theta6.mat <- list()
  Theta6.mat[[d]] <- list()
  for (i in 1:d)
    Theta6.mat[[i]] <- list()
  
  for (i in 1:d)
    for (j in 1:d)
    {  
      temp <- numeric()
      for (k in 1:d)     
        for (ell in 1:d)    
          temp <- rbind(temp, elem(i,d)+2*elem(k,d)+2*elem(ell,d)+elem(j,d))
      
      Theta6.mat[[i]][[j]] <- temp
    }
  
  return(Theta6.mat)
}

default.gridsize <- function(d)
{
  if (d==1)
    gridsize <- 401
  else if (d==2)
    gridsize <- rep(151,d)
  else if (d==3)
    gridsize <- rep(51, d)
  else if (d==4)
    gridsize <- rep(21, d)
  else
    gridsize <- NA
  
  return(gridsize)
}

default.bgridsize <- function(d)
{
  if (d==1)
    gridsize <- 401
  else if (d==2)
    gridsize <- rep(151,d)
  else if (d==3)
    gridsize <- rep(31, d)
  else if (d==4)
    gridsize <- rep(21, d)
  else
    gridsize <- NA
  
  return(gridsize)
}

