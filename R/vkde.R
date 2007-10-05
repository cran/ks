###############################################################################
## Variable bandwidth selectors
###############################################################################


###############################################################################
## Prepare data x for variable KDE
##
## Parameters
## x - data
## type - "abram" - for Abramson
##      - "sain" - for Sain (2002)
##      - "clust" - pre-clustered
##
## Returns
## list of with fields
## x - data
## nclust - number of data points in each cluster 
## in - index of cluster membership
###############################################################################


pre.vkde <- function(x, type, G, pre="scale", small=5, nu, dist="euclidean",
    link="average", stop.rule="dudahart", sig.level=0.001)
{
  typ <- substr(tolower(type),1,1)

  if (typ=="a")
    x.pvk <- pre.abram(x=x)
  else if (typ=="s")
    x.pvk <- pre.sain(x=x, G=G, pre=pre, small=small)
  else if (typ=="c")
    x.pvk <- pre.clust(x=x, nu=nu, dist=dist, link=link, stop.rule=stop.rule,
                      sig.level=sig.level, pre=pre)

  class(x.pvk) <- "pvkde"
  return (x.pvk)
}


###############################################################################
## Plot pre-clustered data 
#
## Parameters
## x.pc - pre-clustered data points
#
## Returns
## Plots data with each cluster denoted 1,2,... and in col=1,2,...
###############################################################################

plot.pvkde <- function(x, ...)
{
  x.pvk <- x
  rm(x)
  plot(x.pvk$x, cex=0, ...)
  text(x.pvk$x[,1], x.vk$x[,2], x.pvk$ind, col=x.pvk$ind)
}



Hlscv.vkde <- function(x, Hstart, pre="scale", univar=FALSE)
{
  if (length(x$nclust)==length(x$x))
    H <- Hlscv.ab(x, pre=pre)
  else
    H <- Hlscv.pc(x, Hstart=Hstart, univar=univar)
}



###############################################################################
### Abramsom selectors 
###############################################################################

###############################################################################
## Prepare data x for Abramson variable KDE 
##############################################################################

pre.abram <- function(x)
{
  if (is.vector(x)) n <- 1
  else n <- nrow(x)
  return(pre.clust(x, nu=n))
}


###############################################################################
## LSCV for Abramson type bandwidth matrices using normal reference pilot KDE
## Parameters
## x - data values
#
## Returns
## H_LSCV,AB
###############################################################################

lscv.mat.ab <- function(x, h, fhatp)
{
  d <- ncol(x)
  n <- nrow(x)
  h.ab <- h^2 * fhatp^(-1) 

  if (is.vector(x))
  {
    n <- 1; d <- 1; x1 <- x[1]; x2 <- x[2]
  }
  else
  {
    n <- nrow(x); x1 <- x[,1]; x2 <- x[,2]
  }
  result <- .Fortran("dmvnorm_2d_sum_ab", as.double(x1), as.double(x2), as.double(h.ab),
               as.double(h.ab), as.integer(n), as.integer(1), as.double(0),
               PACKAGE="ks")
  lscv1 <- result[[7]]

  result <- .Fortran("dmvnorm_2d_sum_ab", as.double(x1), as.double(x2), as.double(h.ab),
               as.double(h.ab), as.integer(n), as.integer(0), as.double(0),
               PACKAGE="ks")
  lscv2 <- result[[7]]
  
  return(lscv1/n^2 - 2/(n*(n-1))*lscv2)
}

###############################################################################
## Find Abramson type bandwidth matrices with LSCV 
## 
## Parameters
## x - data values
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data
#
## Returns
## H_LSCV,A
###############################################################################

Hlscv.ab <- function(x, pre="scale")
{
  d <- 2
  if (substr(pre,1,2)=="sc")
    pre <- "scale"
  else if (substr(pre,1,2)=="sp")
    pre <- "sphere"
  if(!is.matrix(x)) x <- as.matrix(x)

  ## pre-transform data
  if (pre=="sphere")
    x.star <- pre.sphere(x)
  else if (pre=="scale")
    x.star <- pre.scale(x)
  S.star <- var(x.star)
  n <- nrow(x.star)

  G.nr.star <- (4/(n*(d+2)))^(2/(d+4)) * S.star
  fhatp <- kde(x.star, G.nr.star, eval.points=x.star)$est

  lscv.mat.ab.temp <- function(h)
  {
    return(lscv.mat.ab(x.star, h, fhatp))
  }

  hend <- sum(diag((4/(n*(d+2)))^(2/(d+4))*var(x)*5))
  result <- optimise(lscv.mat.ab.temp, interval=c(0, hend))
  H.ab <- result$min^2 * diag(d)

  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
  H.ab <- S12 %*% H.ab %*% S12

  H.lscv.ab <- numeric()
  for (i in 1:nrow(x))
    H.lscv.ab <- rbind(H.lscv.ab, H.ab * fhatp[i]^(-1))
  
  return(H.lscv.ab)
}


###############################################################################
### Pre-clustered selectors
###############################################################################


###############################################################################
## Exact MISE for normal mixtures (for pre-clustered KDE)
#
## Parameters
## mus - means
## Sigmas - variances
## Props - vector of proportions of each mixture component 
## Hpc - bandwidth matrices
## samppc - sample sizes
#
## Returns
## Exact MISE for normal mixtures
###############################################################################

mise.mixt.pc <- function(H.pc, mus, Sigmas, props, samp.pc)
{
  d <- ncol(Sigmas)
  n <- sum(samp.pc)
  m <- length(samp.pc)
  M <- length(props)
  if (M != m)
    stop("Number of bandwidth matrices not equal to number of mixture components") 
 
  mise1 <- 0
  mise2 <- 0
  mise3 <- 0
  zero.mat <- matrix(0, nc=d, nr=d) 

  for (j in 1:m)
  {
    Hj <- H.pc[((j-1)*d+1):(j*d),]
    nj <- samp.pc[j]
    mise1 <- mise1 +  nj* det(Hj)^(-1/2) * (4*pi)^(-d/2)

    for (jprime in 1:m)
    {
      mise21 <- 0 
      if (j != jprime)
      { 
        Hjprime <- H.pc[((jprime-1)*d+1):(jprime*d),]
        mise21 <- mise21 +  Xi(mus, Sigmas, M, Hj, Hjprime)
      }
    }
    mise2 <- mise2 + nj * (nj-1) * mise21
    mise3 <- mise3 + nj * Xi(mus, Sigmas, M, Hj, zero.mat)
  }
  mise4 <- Xi(mus, Sigmas, M, zero.mat, zero.mat)
  
  return(mise1/n^2 + drop(props %*% (mise2/n^2 - 2/n*mise3 + mise4) %*% props))

}

###############################################################################
## ISE for normal mixtures (pre-clustered data)
## 
## Parameters
## x.pc - pre-clustered data values
## H.pc - bandwidth matrices
## mus - matrix of means
## Sigmas - matrix of covariance matrices 
## props - mixing proportions 
#
## Returns
## ISE
###############################################################################

ise.mixt.pc <- function(x.pc, H.pc, mus, Sigmas, props, lower, upper,
                        gridsize, stepsize)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")
  x <- x.pc$x 
  n <- nrow(x)
  m <- length(props)
  d <- ncol(Sigmas)
  ind.lab <- sort(unique(x.pc$ind))
  if (is.vector(mus)) mus <- matrix(mus, nc=d)
  Hs <- numeric(0)
  
  for (i in 1:n)
  {
    clust <- which(x.pc$ind[i]==ind.lab)
    H <- H.pc[((clust-1)*d+1):(clust*d),]
    Hs <- rbind(Hs, H)
  }

  #ise1 <- 0
  #ise2 <- 0
  #ise3 <- 0
  
  #for (i in 1:n)
  #{
  ##  clusti <- which(x.pc$ind[i]==ind.lab)
  ##  Hi <- H.pc[((clusti-1)*d+1):(clusti*d),]
  ##  
  ##  for (ip in 1:n)
  ##  {
  ##    clustip <- which(x.pc$ind[i]==ind.lab)
  ##    Hip <- H.pc[((clustip-1)*d+1):(clustip*d),]
  ##    ise1 <- ise1 + dmvnorm.mixt(x[i,], mus=x[ip,], Sigmas=Hi + Hip, props=1)
  ##  }
   
  ##  for (k in 1:m)
  ##  {
  ##    Sigmak <- Sigmas[((k-1)*d+1):(k*d),]
  ##    ise2 <- ise2 + props[k]*dmvnorm.mixt(x[i,], mus=mus[k,], Sigmas=Hi + Sigmak, props=1)
  ##  }
  #}
  #for (k in 1:m)
  #{
  ##  Sigmak <- Sigmas[((k-1)*d+1):(k*d),]
  ##  for (kp in 1:m)
  ##  {
  ##    Sigmakp <- Sigmas[((kp-1)*d+1):(kp*d),]
  ##    ise3 <- ise3 + props[k]*props[kp]*dmvnorm.mixt(mus[k,], mus[kp,], Sigmakp + Sigmak, props=1)
  ##  }
  #}
  
  #ise <- ise1/n^2 - 2*ise2/n + ise3
  
  if (!missing(gridsize))
  {  
    xx <- seq(lower[1], upper[1], length=gridsize[1])
    yy <- seq(lower[2], upper[2], length=gridsize[2])
  }
  else if (!missing(stepsize))
  {
    xx <- seq(lower[1], upper[1], by=stepsize[1])
    yy <- seq(lower[2], upper[2], by=stepsize[2]) 
  }
 
  xxyy <- permute(list(xx, yy))
  fhat <- dmvnorm.mixt(x=xxyy, mus=x, Sigma=Hs, props=rep(1/n,n))
  mixt <- dmvnorm.mixt(x=xxyy, mu=mus, Sigma=Sigmas, props=props)
  stepsize <- c(xx[1]-xx[2], yy[1]-yy[2])
  
  ise <- sum((fhat-mixt)^2*prod(stepsize))

  return(ise)
}




###############################################################################
## Finds the local mimina in matrix of data values for 2D
#
## Parameters
## heights - data matrix
#
## Returns
## matrix of indices of local minima
###############################################################################


mvpeaks.2d <- function(height, span=3)
{
  ## function by Petr Pikal 
  peaks <- function (series, span = 3)
  {
    z <- embed(series, span)
    s <- span%/%2
    v <- max.col(z) == 1 + s
    result <- c(rep(FALSE, s), v)
    result <- result[1:(length(result) - s)]
    result
  }
  
  nc <- ncol(height)
  nr <- nrow(height)

  ## find modes in x and y directions
  x.max <- vector()
  y.max <- vector()
  for (i in 2:(nc-1))  
    x.max <- cbind(x.max, peaks(height[,i], span=span))
 
  for (i in 2:(nr-1))
    y.max <- rbind(y.max, peaks(height[i,], span=span))

  ## combine modes for each direction
  modes.ind <- which(x.max & y.max, arr.ind=TRUE)

  return(modes.ind)
}


###############################################################################
## Finds the nearest mode to data point
#
## Parameters
## x - data matrix
## modes - matrix of sample modes
#
## Returns
## vector of indices of nearest mode 
###############################################################################

nearest.mode <- function(x, modes)
{
  dist <- matrix(0, nr=nrow(x), nc=nrow(modes))
  lab <- rep(0, length=nrow(x))
  for (i in 1:nrow(x))
  {
    for (j in 1:nrow(modes))
      dist[i,j] <- sqrt(sum((x[i,] - modes[j,])^2))
    lab[i] <- which.min(dist[i,])
  }
  return(lab)
}

###############################################################################
## Prepare data points for Sain's partitioned b/w selector  
#
## Parameters
## x - data matrix
## G - pilot b/w matrix
## pre - pre-transformation
## small - smallest number of points allowed to form a separate cluster
#
## Returns
## list of values with components (suitable for use with Hlscv.pc)
## x - data values
## ind - vector of indices indicating which partition class
## nclust - vector of number of points in each partition class 
###############################################################################

pre.sain <- function(x, G, pre="scale", small=5)
{
  d <- 2
  if (substr(pre,1,2)=="sc")
    pre <- "scale"
  else if (substr(pre,1,2)=="sp")
    pre <- "sphere"
  if(!is.matrix(x)) x <- as.matrix(x)

  ## pre-transform data
  if (pre=="sphere")
    x.star <- pre.sphere(x)
  else if (pre=="scale")
    x.star <- pre.scale(x)
  S.star <- var(x.star)
  n <- nrow(x.star)

  ## pilot KDE
  if (missing(G)) G <- (4/(n*(d+2)))^(2/(d+4)) * S.star
  fhatp <- kde(x.star, G, gridsize=c(100,100))

  ## find indices for sample modes for pilot KDE
  modes.ind <- mvpeaks.2d(fhatp$estimate)
 
  ## find value of modes from indices
  if (is.vector(modes.ind))
    modes <- c(fhatp$eval.points[[1]][modes.ind[1]],
               fhatp$eval.points[[2]][modes.ind[2]])
  else
  {
    modes <- numeric()
    for (i in 1:nrow(modes.ind))
       modes <- rbind(modes, c(fhatp$eval.points[[1]][modes.ind[i,1]],
                             fhatp$eval.points[[2]][modes.ind[i,2]]))
  }

  ## associate data points to nearest mode
  if (is.vector(modes))
  { 
    modes.lab <- rep(1, length=nrow(x))
    nclust <- nrow(x)
  }
  else
  {
    modes.lab <- nearest.mode(x.star, modes)
    mlab <- sort(unique(modes.lab))
    nclust <- rep(0, length(mlab)) 
    for (i in 1:length(mlab))
      nclust[i] <- sum(modes.lab==mlab[i])
  }

  if (sum(nclust<=small)>0)
  {
    small.ind <- which(nclust<=small)
    modes.ind <- modes.ind[-small.ind,]
    
    if (is.vector(modes.ind))
      modes <- c(fhatp$eval.points[[1]][modes.ind[1]],
               fhatp$eval.points[[2]][modes.ind[2]])
    else
    {
      modes <- numeric()
      for (i in 1:nrow(modes.ind))
        modes <- rbind(modes, c(fhatp$eval.points[[1]][modes.ind[i,1]],
                             fhatp$eval.points[[2]][modes.ind[i,2]]))
    }
    if (is.vector(modes))
    { 
      modes.lab <- rep(1, length=nrow(x))
      nclust <- nrow(x)
    }
    else
    {
      modes.lab <- nearest.mode(x.star, modes)
      mlab <- sort(unique(modes.lab))
      nclust <- rep(0, length(mlab)) 
      for (i in 1:length(mlab))
        nclust[i] <- sum(modes.lab==mlab[i])
    }
  }
  #if (plot.it)
  #{
  ##  plot.kde(fhatp, cont=seq(0,100, by=10), col="grey", cex=0)
  ##  text(x.star[,1], x.star[,2], modes.lab)
  ##  if (is.vector(modes))
  ##    text(modes[1], modes[2], col=2)
  ##  else
  ##    text(modes, col=2)
  #}
  
  return(list(x=x, ind=modes.lab, nclust=nclust))
}


###############################################################################
## Pre-cluster data set into nu clusters
#
## Parameters
## x - data points
## nu - number of clusters
## dist - distance function
## link - linking function
## stop.rule  - stopping rule
## sig.level - (step-wise only) significance level
#
## Returns
## List with elements
## x - data matrix
## ind - index labels to indicate cluster membership
## nclust - number of members in each cluster
###############################################################################

pre.clust <- function(x, nu, dist="euclidean", link="average",
                      stop.rule="dudahart", sig.level=0.001, pre="scale")
{
  n <- nrow(x)
  
  ## pre-transform data before clustering
  if (pre=="sphere")
    xstar <- pre.sphere(x)
  else
    xstar <- pre.scale(x)
  
  xd <- dist(xstar, method=dist)
  dendo <- hclust(xd, method=link)

  ## find optimal number of clusters according to stopping rule
  if (missing(nu))
    nu <- optnum.clust(xstar, dendo, stop.rule=stop.rule, sig.level=sig.level)
  clust.ind <- cutree(dendo, k=nu)
  nclust <- rep(0, nu) 

  for (i in 1:nu)
    nclust[i] <- length(which(clust.ind==i))

  return(list(x=x, ind=clust.ind, nclust=nclust))
}



###############################################################################
## Sum of squares decomposition
#
## Parameters
## x.pc - data matrix
#
## Returns
## List with elements
## B - between clusters sum of squares
## W - within clusters sum of squares
###############################################################################

ssq.decomp <- function(x.pc)
{
  x <- x.pc$x
  nu <- length(x.pc$nclust)
  d <- ncol(x.pc$x)
  
  xbar <- apply(x, 2, mean)
  xbari <- matrix(0, nc=d, nrow=nu)
  clust.ind <- x.pc$ind
  nclust <- x.pc$nclust 
  clust.indg <- sort(unique(clust.ind))

  ssq <- function(x, w=rep(1,nrow(x)))
  {
    ssq <- 0
    for (i in 1:nrow(x))
      ssq <- ssq + w[i] *(x[i,] %*% t(x[i,]))
    
    return(ssq)
  }
  
  xW <- numeric()   ## xW is x (in cluster) - cluster mean 
  xB <- numeric()   ## XB is x (in cluster) - total mean
  count <- 0

  for (i in clust.indg)
  {
    count <- count + 1
    xpci <- x[clust.ind==i,]

    ## more than one element in cluster
    if (is.matrix(xpci))
    {
      xbari <- apply(xpci, 2, mean)
      xB <- rbind(xB, xbari - xbar)
      for (j in 1:nclust[count])
        xW <- rbind(xW, xpci[j,] - xbari)
    }
    ## only one element in cluster
    else
    {
      xW <- rbind(xW, rep(0, d))
      xB <- rbind(xB, xpci - xbar)
    }
  }
  W <- ssq(xW)
  B <- ssq(xB, x.pc$nclust)
 
  return(list(B=B,W=W))
}  
 

###############################################################################
## Find optimal number of clusters in data set
#
## Parameters
## x - data points
## dendo - dendogram (output from hclust)
## stop.rule - stopping rule
## sig.level - (step-wise only) significance level of test
#
## Returns
## Optimal number of clusters
###############################################################################

optnum.clust <- function(x, dendo, stop.rule, sig.level)
{
  d <- ncol(x)
  n <- nrow(x)
  tree.temp <- clust.tree(dendo)
  tree <- tree.temp$tree
  ## tree contains all info needed to reconstruct hierarchy
  tree.ind <- tree.temp$ind
  nu <- 1 ## nu is number of clusters 

  ## initial test - compare SS for 1 cluster with that for 2 clusters
  i <-  nrow(tree)
  W1 <- var(x) * (n-1)
  tree.ind2 <- sort(tree[i,1:2])
  x.ind2 <- (tree.ind[,i]==tree.ind2[1]) | ( tree.ind[,i]==tree.ind2[2])
  wti1 <- which(tree.ind[,i]==tree.ind2[1])
  wti2 <- which(tree.ind[,i]==tree.ind2[2])
  nc2 <- c(length(wti1), length(wti2))
  x.pc2 <- list(x=x[sort(c(wti1, wti2)),], ind=tree.ind[x.ind2,i], nclust=nc2)
  W2 <- ssq.decomp(x.pc2)$W
  trW1 <- sum(diag(W1))
  trW2 <- sum(diag(W2))

  if (stop.rule=="dudahart")
  {
    teststat <- trW2/trW1
    nc1 <- n
    critval <- 1-2/(pi*d) - qnorm(1-sig.level)*sqrt(2*(1 - 8/(pi^2*d))/(nc1*d))
 
    while ((teststat < critval) & (nu < n))
    {
      i <- i - 1   ## move backwards up tree
      nu <- nu + 1 ## increment number of clusters 
      ind1 <- tree[i,3]
      tree.ind2 <- sort(tree[i,1:2])
      nc1 <- length(which(tree.ind[,i+1]==ind1)) ## number in cluster 1
      wti1 <- which(tree.ind[,i]==tree.ind2[1])  ## indices of elem in clust 1 
      wti2 <- which(tree.ind[,i]==tree.ind2[2])  ## indices of elem in clust 2
      x.ind1 <- tree.ind[,i+1]==ind1

      ## reconstruct cluster 1 x.pc1 and cluster 2 x.pc2 
      x.pc1 <- list(x=x[which(x.ind1),], ind=tree.ind[x.ind1,i+1], nclust=nc1)
      x.ind2 <- (tree.ind[,i] ==tree.ind2[1]) | (tree.ind[,i]==tree.ind2[2])
      nc2 <- c(length(wti1), length(wti2))
      x.pc2 <- list(x=x[sort(c(wti1, wti2)),], ind=tree.ind[x.ind2,i],
                    nclust=nc2)
      
      W1 <- ssq.decomp(x.pc1)$W
      W2 <- ssq.decomp(x.pc2)$W
      trW1 <- sum(diag(W1))
      trW2 <- sum(diag(W2))
      teststat <- trW2/trW1
      critval <- 1-2/(pi*d) - qnorm(1-sig.level)*sqrt(2*(1-8/(pi^2*d))/(nc1*d))
    }
  }

  return(nu)
}

###############################################################################
## Reconstructs clustering hierarchy 
#
## Parameters
## dendo - output from hclust (of n data points)
#
## Returns
#
## (n-1) x 3 matrix with i-th row: 
## Columns 1, 2 - index label (output from cutree(dendo, k=i)) of the two
##                clusters being merged
## Column 3 - index label of merged cluster (output from cutree(dendo, k=i+1))
###############################################################################

clust.tree <- function(dendo)
{ 
  merge <- dendo$merge
  ## merge is what we need to reconstruct hierarchy
  n <- nrow(merge) + 1
  clust.ind <- as.matrix(cutree(dendo, k=n:1))
  clust.indn <- clust.ind[,-n]
  index <- matrix(0, nrow=n-1, ncol=3)

  ## modify 'merge' to format which 'optnum.clust' requires as input 
  for (i in 1:(n-1))
  {
    m <- merge[i,]
    if (m[1] < 0)
      index[i,1] <- clust.indn[-m[1],i]
    else if (m[1] > 0)
      index[i,1] <- clust.indn[which(clust.indn[,m[1]]==index[m[1],3]),i][1]
    
    if (m[2] < 0)
      index[i,2] <- clust.indn[-m[2],i]
    else if (m[2] > 0)
      index[i,2] <- clust.indn[which(clust.indn[,m[2]]==index[m[2],3]),i][1]
    
    index[i,3] <- min(index[i,1:2]) 
  }

  return(list(tree=index, ind=clust.ind))
}


###############################################################################
## Pre-clustered bandwidth selectors
#
## Parameters
## x.pc - pre-clustered data points
## H.pc - bandwidth matrices
#
## Returns
## Bandwidth matrix for each cluster
###############################################################################

lscv.mat.pc <- function(x.pc, H.pc)
{
  n <- sum(x.pc$nclust)
  M <- length(x.pc$nclust)
  d <- ncol(x.pc$x)
  H <- list()
  ind.lab <- sort(unique(x.pc$ind))
  for (j in 1:M)
  {  
    H[[j]] <- H.pc[((j-1)*d+1):(j*d),]
    while (det(H[[j]]) < 1e-12)
      for (i in 1:d) H[[j]][i,i] <- 1.01* H[[j]][i,i]
  }
    
  lscv1 <- 0
  lscv2 <- 0
  
  for (j in 1:M)
  {
    xj <- x.pc$x[x.pc$ind==ind.lab[j],]
    for (jprime in 1:M)
    {
      xjprime <- x.pc$x[x.pc$ind==ind.lab[jprime],]
      lscv1 <- lscv1 + dmvnorm.2d.sum.pc(xj,xjprime, H[[j]]+H[[jprime]],inc=1)
      lscv2 <- lscv2 + dmvnorm.2d.sum.pc(xj,xjprime, H[[jprime]], inc=1)
    }
  }
  e <- 0
  for (i in 1:n)
    e <- e + dmvnorm(rep(0,d),  rep(0,d), H[[which(ind.lab==x.pc$ind[i])]])
  lscv2 <- lscv2 - e

  return(lscv1/n^2 - 2/(n*(n-1))*lscv2)     
}


###############################################################################
## Finds the bandwidth matrices that minimises LSCV 
## 
## Parameters
## x.pc - pre-clustered data values
## Hstart - initial bandwidth matrices
#
## Returns
## H_PC,LSCV
###############################################################################

Hlscv.pc <- function(x.pc, Hstart, univar=FALSE)
{
  x <- x.pc$x
  n <- nrow(x)
  d <- ncol(x)
  dstar <- d*(d+1)
  nu <- length(x.pc$nclust)
  RK <- (4*pi)^(-d/2)
  ind.lab <- sort(unique(x.pc$ind))
  
  ## initial conditions based on normal reference b/w within each cluster
  if (missing(Hstart))
  {
    Hstart <- numeric()
    for (j in 1:nu)
    {
      xj <- x[x.pc$ind==ind.lab[j],]

      if (!univar)
        if (is.vector(xj))
          Hstart <- rbind(Hstart,sqrt(4/(n*(d + 2))^(2/(d + 4)))*diag(d))
        else
          Hstart <- rbind(Hstart,matrix.sqrt((4/(n*(d + 2)))^(2/(d + 4))*var(xj)))
      else  
        #if (is.vector(xj))
          Hstart <- c(Hstart,sqrt(4/(n*(d + 2))^(2/(d + 4))))
        #else 
        ##  Hstart <- c(Hstart, sqrt((4/(n*(d + 2)))^(2/(d + 4))*max(abs(var(xj)))))
    }                 
  }

  if (univar)
  {
    x.pc2 <- x.pc
    x.pc2$x <- pre.scale(x.pc$x)
    lscv.mat.pc.temp <- function(hpc)
    {
      H.pc <- numeric()
      for (j in 1:nu)
        H.pc  <- rbind(H.pc,  hpc[j]*hpc[j]*diag(d))
        ##  ensures that H.pc are positive definite
       
      return(lscv.mat.pc(x.pc2, H.pc))
    }
    if (nu==1)
    {
      hend <- sqrt(sum(diag((4/(n*(d+2)))^(2/(d+4))*var(x))))
      result <- optimise(lscv.mat.pc.temp, interval=c(0,10*hend))
      H.pc <- result$min^2 * diag(d)* diag(sd(x)^2)
    }
    else
    {  
      result <- optim(Hstart, lscv.mat.pc.temp, method="Nelder-Mead")
      H.pc <-  matrix(0, nc=d, nr=nu*d) 
      for (j in 1:nu)
        H.pc[((j-1)*d+1):(j*d),] <- result$par[j]^2 * diag(d) * diag(sd(x)^2)
    }
  }
  else    
  {
    lscv.mat.pc.temp <- function(vechHpc)
    {
      H.pc <- numeric()
      for (j in 1:nu)
      {   
        vechHj <- vechHpc[((j-1)*1/2*dstar+1): (j*1/2*dstar)]
        H.pc  <- rbind(H.pc, invvech(vechHj) %*% invvech(vechHj))
        ##  ensures that H.pc are positive definite
      }
      
      return(lscv.mat.pc(x.pc, H.pc))
    }
    result <- optim(vech.cat(Hstart), lscv.mat.pc.temp, method="Nelder-Mead")
    H.pc <- invvech.cat(result$par,d)
    for (j in 1:nu)
      H.pc[((j-1)*d+1):(j*d),] <- H.pc[((j-1)*d+1):(j*d),] %*% H.pc[((j-1)*d+1):(j*d),]
  }
 
  return(H.pc)
}



