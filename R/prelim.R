
#############################################################################
## Basic vectors and matrices and their operations
#############################################################################
  
## Vec operator 

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

## Vech operator

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

  
## Inverse vec operator 
    
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

## Inverse vech operator 

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

## Trace of matrix

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


## Elementary vector 
    
elem <- function(i, d)
{
  elem.vec <- rep(0, d)
  elem.vec[i] <- 1
  
  return(elem.vec)
}      

## Commutation matrix (taken from MCMCglmmm library)

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


differences <- function(x, y, upper=FALSE, ff=FALSE, Kpow=0)
{
  if (missing(y)) y <- x
  if (is.vector(x)) x <- t(as.matrix(x))
  if (is.vector(y)) y <- t(as.matrix(y))

  nx <- nrow(x)
  ny <- nrow(y)
  d <- ncol(x)

  if (ff) difs <- ff(init=0, dim=c(nx*ny,d))
  else difs <- matrix(ncol=d,nrow=nx*ny)

  for (j in 1:d)
  {
    difs[,j] <- rep(x[,j], times=ny) - rep(y[1:ny,j], each=nx)
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

DrL0<-function(d,r,Sdr.mat, Sdr.flag=TRUE, verbose=FALSE, thin=1)
{
  if (!Sdr.flag)
  {
    v <- Kpow(A=vec(diag(d)),pow=r/2)
    DL0<-(-1)^(r/2)*(2*pi)^(-d/2)*OF(r)*matrix(Sdrv(d=d, r=r, v=v, thin=thin, verbose=verbose), ncol=1) 
  }
  else 
  { 
    if (missing(Sdr.mat)) Sdr.mat <- Sdr(d=d, r=r)
    DL0<-(-1)^(r/2)*(2*pi)^(-d/2)*OF(r)*(Sdr.mat%*%Kpow(A=vec(diag(d)),pow=r/2))
  }
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




#### Permutations with repetitions of the first d naturals (1:d) taking k elements at a time
#### There are d^k of them, each having length k => We arrange them into a matrix of order d^k times k
#### Each row represents one permutation
#### Second version: filling in the matrix comlumn-wise (slightly faster)
    
perm.rep <-function(d,k)
{
    PM <- matrix(nrow=d^k,ncol=k)
    for(pow in 0:(k-1)){
        t2<-d^pow
        p1<-1
        while(p1<=d^k){
            for(al in 1:d){for(p2 in 1:t2){
                PM[p1,k-pow]<-al
                p1<-p1+1}}}}
    return(PM)
}


### Symmetriser applied to a vector

Sdrv <- function(d,r,v, thin=1, verbose=FALSE)
{   
    if(length(v)!=d^r){stop("The length of v must equal d^r")}
    per.rep<-perm.rep(d,r)
    nper<-factorial(r)
    ##nper.rep<-d^r
    dpow<-d^((r-1):0)
    dpow.mat<-matrix(rep(dpow,d^r),ncol=r,nrow=d^r,byrow=TRUE)    
    if (verbose) pb <- txtProgressBar()
    
    ## modified from the permn function in the combinat library 02/2012   
  
    Sv <- rep(0,d^r)
    x <- seq(r)
    n <- length(x)
    ##nofun <- is.null(fun)
    ##out <- vector("list", gamma(n + 1))
    p <- ip <- seqn <- 1:n
    d <- rep(-1, n)
    d[1] <- 0
    m <- n + 1
    p <- c(m, p, m)
    i <- 1
    nper.thin <- 0
    use <- -c(1, n + 2)
   
    while (m != 1) {
      ##out <- if (nofun) x[p[use]] else fun(x[p[use]], ...)
      out <- x[p[use]]
      if (i%%thin == 0 | i==1)
      { 
        ## next 2 lines added to compute Sdr%*%v product
        if (verbose) setTxtProgressBar(pb, i/nper)
        ##cat(paste(i,"/", nper, "\n", sep=""))
        Sv <- Sv + v[1+rowSums((per.rep[,out]-1)*dpow.mat)]
        nper.thin <- nper.thin+1
      }

      i <- i + 1
      m <- n
      chk <- (p[ip + d + 1] > seqn)
      m <- max(seqn[!chk])
      if (m < n) d[(m + 1):n] <- -d[(m + 1):n]
      index1 <- ip[m] + 1
      index2 <- p[index1] <- p[index1 + d[m]]
      p[index1 + d[m]] <- m
      tmp <- ip[index2]
      ip[index2] <- ip[m]
      ip[m] <- tmp
   }
   permnr <- Sv/nper.thin  
   if (verbose) close(pb)
   return(permnr)

    ## attempt at speeding up the loop over permutations
    ##per.indices <- block.indices(nrow(per.rep)*ncol(per), nrow(per), block.limit=1e7)
    ##Sv <- 0
    ##for (i in 1:(length(per.indices)-1))
    ##{ 
    ##  if (verbose) setTxtProgressBar(pb, i/(length(per.indices)-1))
    ##  nper.temp <- diff(per.indices[i:(i+1)])
    ##  Sv.mat <- (per.rep[,t(per.mat[per.indices[i]:(per.indices[(i+1)]-1),])]-1)*matrix(rep(dpow,d^r),ncol=r*nper.temp,nrow=d^r,byrow=TRUE) 
    ##  Sv <- Sv + rowSums(matrix(v[1+apply(array(Sv.mat, dim=c(d^r,r,nper.temp)), 3, rowSums)], nrow=d^r))
    ##}
}

## default block indices for double sums
block.indices <- function(nx, ny, d, r=0, diff=FALSE, block.limit=1e6, npergroup)
{
  if (missing(npergroup)) 
  { 
    if (diff) npergroup <- max(c(block.limit %/% (nx*d^r), 1))
    else npergroup <- max(c(block.limit %/% nx,1))
  }
  nseq <- seq(1, ny, by=npergroup)
  if (tail(nseq,n=1) <= ny) nseq <- c(nseq, ny+1)
  if (length(nseq)==1) nseq <- c(1, ny+1)
  return(nseq)
}
