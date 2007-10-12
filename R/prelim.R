
###############################################################################
# Basic operators and functions
###############################################################################
    
###############################################################################
# Vec operator 
#
# Takes square matrix x and returns its elements stacked column wise in a
# vector
###############################################################################

vec <- function(x)
{
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
  d <- ncol(x)
  vechx <- vector()
  
  for (j in 1:d)
    vechx <- c(vechx, x[j:d,j])
                   
  return(vechx)
}

# Analogue for stacked matrix
vech.cat <- function(x)
{
  d <- ncol(x)
  num <- nrow(x)/ncol(x)
  vechx <- vector()
  
  for (j in 1:num)
    vechx <- c(vechx, vech(x[((j-1)*d+1) : (j*d),]))
                   
  return(vechx)
}

  
###############################################################################
# Inverse vec operator 
#
# Takes vector x and stacks its elements into columns of a matrix. 
###############################################################################
    
invvec <- function(x)
{
  d <- sqrt(length(x))
  invvecx <- matrix(0, nr = d, nc = d)
  
  for (j in 1:d)
    invvecx[,j] <- x[c(1:d) + (j-1)*d]
  
  return(invvecx)
}

###############################################################################
# Inverse vech operator 
#
# Takes vector x and stacks its elements into a symmetric matrix 
###############################################################################

invvech <- function(x)
{
  d <- (-1 + sqrt(8*length(x) + 1))/2
  invvechx <- matrix(0, nr = d, nc = d)

  for (j in 1:d)
    invvechx[j:d,j] <- x[1:(d-j+1)+ (j-1)*(d - 1/2*(j-2))]
  
  invvechx <- invvechx + t(invvechx) - diag(diag(invvechx))
  return(invvechx)
}

invvech.cat <- function(x, d)
{
  dstar <- d*(d+1)/2
  num <- length(x)/dstar
  invvechx <- numeric()
  for (j in 1:num)
    invvechx <- rbind(invvechx, invvech(x[((j-1)*dstar+1) : (j*dstar)]))

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
  orep <- final.len <- prod(sapply(args, length))
  for (i in 1:nargs) {
    x <- args[[i]]
    nx <- length(x)
    orep <- orep/nx
    cargs[[i]] <- rep(rep(x, rep(rep.fac, nx)), orep)
    rep.fac <- rep.fac * nx
  }
  do.call("cbind", cargs)
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

###############################################################################
# Matrix square root - taken from Stephen Lake 
# http://www5.biostat.wustl.edu/s-news/s-news-archive/200109/msg00067.html
###############################################################################

matrix.sqrt <- function(A)
{
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

###############################################################################
# For a list of matrices, this returns the dimensions of each matrix
###############################################################################

list.length <- function(x)
{
  ell <- length(x)
  len <- matrix(0, nr=ell, nc=2)
  for (i in 1:ell)
    len[i,] <- dim(x[[i]])

  return(len)
}

###############################################################################
# Pre-sphering
# Parameters
# x - data points
#
# Returns
# Pre-sphered x values
###############################################################################

pre.sphere <- function(x)
{
  S <- var(x)
  Sinv12 <- matrix.sqrt(chol2inv(chol(S)))
  
  x.sphered <- matrix(0, nc=ncol(x), nr=nrow(x))
  for (i in 1:nrow(x))
    x.sphered[i,] <- Sinv12 %*% x[i,]    

  return (x.sphered)
}

pre.sphere.pc <- function(x.pc)
{
  g <- length(x.pc$nclust)
  d <- ncol(x.pc$x)
  x.pc1 <- x.pc
  x.pc1$x <- matrix(0, nc=d, nr=nrow(x.pc$x))
  for (j in 1:g)
  {
    xj <- x.pc$x[x.pc$ind==j,]
    x.pc1$x[x.pc$ind==j,] <- pre.sphere(xj) 
  }
 
  return (x.pc1)          
}  

###############################################################################
# Pre-scaling
# Parameters
# x - data points
#
# Returns
# Pre-scaled x values
###############################################################################

pre.scale <- function(x)
{
  x.scaled <- numeric()
  x.sd <- apply(x, 2, sd)
  d <- ncol(x)

  for (i in 1:d)
    x.scaled <- cbind(x.scaled, x[,i] /x.sd[i])
    

  return (x.scaled)
}

pre.scale.pc <- function(x.pc)
{
  g <- length(x.pc$nclust)
  d <- ncol(x.pc$x)
  x.pc1 <- x.pc
  x.pc1$x <- matrix(0, nc=d, nr=nrow(x.pc$x))
  for (j in 1:g)
  {
    xj <- x.pc$x[x.pc$ind==j,]
    x.pc1$x[x.pc$ind==j,] <- pre.scale(xj) 
  }
 
  return (x.pc1)          
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

  for (i in 1:nrow(mat))
    if (identical(r, mat[i,])) return(i)

  return(NA)  
}

is.even <- function(x)
{
  y <- x[x>0] %%2
  return(identical(y, rep(0, length(y))))
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






