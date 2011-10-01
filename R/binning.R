
### Default grid sizes

default.gridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(51, d)
  else if (d==4) gridsize <- rep(21, d)
  else gridsize <- NA
  
  return(gridsize)
}

default.bgridsize <- function(d)
{
  if (d==1)      gridsize <- 401
  else if (d==2) gridsize <- rep(151,d)
  else if (d==3) gridsize <- rep(31, d)
  else if (d==4) gridsize <- rep(21, d)
  else  gridsize <- NA
  
  return(gridsize)
}


########################################################################
## Linear binning
## Courtesy of M Wand 2005
## Extended by T Duong to 3- and 4-dim 2006
########################################################################



binning <- function(x, H, h, bgridsize, xmin, xmax, supp=3.7, w)
{
  x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  if (missing(w)) w <- rep(1,n)
  if (d>1 & !missing(H))
     if (!identical(diag(diag(H)), H))
       stop("Binning requires diagonal bandwidth matrix")
  
  if (missing(h)) h <- rep(0,d)
  if (!missing(H)) h <- sqrt(diag(H))

  if (missing(bgridsize)) bgridsize <- default.gridsize(d)
  if (!(missing(xmin) & missing(xmax)))
  {
    range.x <- list()
    for (i in 1:d)
      range.x[[i]] <- c(xmin[i], xmax[i])
  }
  else
  {
    range.x <- list()
    for (i in 1:d)
      range.x[[i]] <- c(min(x[,i]) - supp*h[i], max(x[,i]) + supp*h[i])
  }
  a <- unlist(lapply(range.x,min))
  b <- unlist(lapply(range.x,max))

  gpoints <- list()
  for (id in 1:d)
    gpoints[[id]] <- seq(a[id],b[id],length=bgridsize[id])  
 
  if (d==1) counts <- linbin.ks(x,gpoints[[1]], w=w) 
  if (d==2) counts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]], w=w)
  if (d==3) counts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]], w=w)
  if (d==4) counts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]], w=w)
 
  bin.counts <- list(counts=counts, eval.points=gpoints, w=w)
  if (d==1) bin.counts$eval.points <- gpoints[[1]]
     
  return(bin.counts)
}



######### R-function:dfltCounts  #########
 
# Obtain default set of grid counts from a 
# multivariate point cloud 'x'.

# Last changed: 18 JUL 2005

dfltCounts <- function(x,gridsize=rep(64,NCOL(x)),h=rep(0,NCOL(x)), supp=3.7, range.x, w)
{
   x <- as.matrix(x)
   d <- ncol(x)
   n <- nrow(x)
   if (missing(w)) w <- rep(1,n)
   
   if (missing(range.x))
   {
     range.x <- list()
     for (id in 1:d)
       range.x[[id]] <- c(min(x[,id])-supp*h[id],max(x[,id])+supp*h[id])  
   }

   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))

   gpoints <- list()
   for (id in 1:d)
      gpoints[[id]] <- seq(a[id],b[id],length=gridsize[id])  
 
   if ((d!=1)&(d!=2)&(d!=3)&(d!=4)) stop("currently only for d=1,2,3,4")

   if (d==1) gcounts <- linbin.ks(x,gpoints[[1]], w=w) 
   if (d==2) gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]], w=w)
   if (d==3) gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]], w=w)
   if (d==4) gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]], w=w)
   
   return(list(counts=gcounts,range.x=range.x))
}

######## End of dfltCounts ########


########################################################################
## Linear binning
########################################################################

linbin.ks <- function(X, gpoints, w)
{
   n <- length(X)
   M <- length(gpoints)
   if (missing(w)) w <- rep(1, n)
   a <- gpoints[1]
   b <- gpoints[M]
   xi <- .C("massdist", x=as.double(X), xmass=as.double(w), nx=as.integer(n),
            xlo=as.double(a), xhi=as.double(b), y=double(M), ny=as.integer(M),
            PACKAGE="stats")$y
   return(xi)
}


linbin2D.ks <- function(X, gpoints1, gpoints2, w)
{
   n <- nrow(X)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist2d", x1=as.double(X[,1]), x2=as.double(X[,2]), n=as.integer(n),
             a1=as.double(a1), a2=as.double(a2), b1=as.double(b1), b2=as.double(b2),
             M1=as.integer(M1), M2=as.integer(M2), weight=as.double(w), est=double(M1*M2), PACKAGE="ks")
   xi <- matrix(out$est, nrow=M1, ncol=M2)

   ## adjust binning weights for boundary points
   xbmaxpos1 <- which(X[,1]>=max(gpoints1))
   xbmaxpos2 <- which(X[,2]>=max(gpoints2))
   if (length(xbmaxpos1)>0)
   {
     ind1 <- findInterval(X[xbmaxpos1,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos1,2], gpoints2)
     xi[ind1, ind2] <- xi[ind1, ind2] + w[xbmaxpos1]
   }
   if (length(xbmaxpos2)>0)
   {
     ind1 <- findInterval(X[xbmaxpos2,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos2,2], gpoints2)
     xi[ind1, ind2] <- xi[ind1, ind2] + w[xbmaxpos2]
   }
   
   return(xi)
}

linbin3D.ks <- function(X, gpoints1, gpoints2, gpoints3, w)
{
   n <- nrow(X)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist3d", x1=as.double(X[,1]), x2=as.double(X[,2]), x3=as.double(X[,3]), n=as.integer(n),
             a1=as.double(a1), a2=as.double(a2), a3=as.double(a3), b1=as.double(b1), b2=as.double(b2), b3=as.double(b3),
             M1=as.integer(M1), M2=as.integer(M2), M3=as.integer(M3), weight=as.double(w), est=double(M1*M2*M3), PACKAGE="ks")
   xi <- array(out$est, dim=c(M1,M2,M3))
   
   ## adjust binning weights for boundary points
   xbmaxpos1 <- which(X[,1]>=max(gpoints1))
   xbmaxpos2 <- which(X[,2]>=max(gpoints2))
   xbmaxpos3 <- which(X[,3]>=max(gpoints3))

   if (length(xbmaxpos1)>0)
   {
     ind1 <- findInterval(X[xbmaxpos1,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos1,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos1,3], gpoints3)
     xi[ind1,ind2,ind3] <- xi[ind1,ind2,ind3] + w[xbmaxpos1]
   }
   if (length(xbmaxpos2)>0)
   {
     ind1 <- findInterval(X[xbmaxpos2,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos2,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos2,3], gpoints3)
     xi[ind1,ind2,ind3] <- xi[ind1,ind2,ind3] + w[xbmaxpos2]
   }
   if (length(xbmaxpos3)>0)
   {
     ind1 <- findInterval(X[xbmaxpos3,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos3,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos3,3], gpoints3)
     xi[ind1,ind2,ind3] <- xi[ind1,ind2,ind3] + w[xbmaxpos3]
   }
   return(xi)
}

linbin4D.ks <- function(X, gpoints1, gpoints2, gpoints3, gpoints4, w)
{
   n <- nrow(X)
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   M4 <- length(gpoints4)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   a4 <- gpoints4[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   b4 <- gpoints4[M4]
   if (missing(w)) w <- rep(1, n)

   ## binning for interior points
   out <- .C("massdist4d", x1=as.double(X[,1]), x2=as.double(X[,2]), x3=as.double(X[,3]), x4=as.double(X[,4]), n=as.integer(n),
             a1=as.double(a1), a2=as.double(a2), a3=as.double(a3), a4=as.double(a4),
             b1=as.double(b1), b2=as.double(b2), b3=as.double(b3), b4=as.double(b4),
             M1=as.integer(M1), M2=as.integer(M2), M3=as.integer(M3), M4=as.integer(M4),
             weight=as.double(w), est=double(M1*M2*M3*M4), PACKAGE="ks")
   xi <- array(out$est, dim=c(M1,M2,M3,M4))

   ## adjust binning weights for boundary points
   xbmaxpos1 <- which(X[,1]>=max(gpoints1))
   xbmaxpos2 <- which(X[,2]>=max(gpoints2))
   xbmaxpos3 <- which(X[,3]>=max(gpoints3))
   xbmaxpos4 <- which(X[,4]>=max(gpoints4))
   
   if (length(xbmaxpos1)>0)
   {
     ind1 <- findInterval(X[xbmaxpos1,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos1,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos1,3], gpoints3)
     ind4 <- findInterval(X[xbmaxpos1,4], gpoints4)
     xi[ind1,ind2,ind3,ind4] <- xi[ind1,ind2,ind3,ind4] + w[xbmaxpos1]
   }
   if (length(xbmaxpos2)>0)
   {
     ind1 <- findInterval(X[xbmaxpos2,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos2,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos2,3], gpoints3)
     ind4 <- findInterval(X[xbmaxpos2,4], gpoints4)
     xi[ind1,ind2,ind3,ind4] <- xi[ind1,ind2,ind3,ind4] + w[xbmaxpos2]
   }
   if (length(xbmaxpos3)>0)
   {
     ind1 <- findInterval(X[xbmaxpos3,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos3,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos3,3], gpoints3)
     ind4 <- findInterval(X[xbmaxpos3,4], gpoints4)
     xi[ind1,ind2,ind3,ind4] <- xi[ind1,ind2,ind3,ind4] + w[xbmaxpos3]
   }
   if (length(xbmaxpos4)>0)
   {
     ind1 <- findInterval(X[xbmaxpos4,1], gpoints1)
     ind2 <- findInterval(X[xbmaxpos4,2], gpoints2)
     ind3 <- findInterval(X[xbmaxpos4,3], gpoints3)
     ind4 <- findInterval(X[xbmaxpos4,4], gpoints4)
     xi[ind1,ind2,ind3,ind4] <- xi[ind1,ind2,ind3,ind4] + w[xbmaxpos4]
   }
   
   return(xi)
}


### lingbin*.ksf are based on Fortrtan code of M. Wand & T. Duong 

linbin.ksf <- function(X, gpoints, truncate=FALSE, w)
{
   n <- length(X)
   M <- length(gpoints)  
   trun <- 0
   if (truncate) trun <- 1
   if (missing(w)) w <- rep(1, n)
     
   a <- gpoints[1]
   b <- gpoints[M]
   out <- .Fortran("linbin",as.double(X),as.integer(n),
           as.double(a),as.double(b),as.integer(M),
           as.integer(trun), as.double(w), double(M),PACKAGE="ks")
   return(out[[8]])
}

linbin2D.ksf <- function(X, gpoints1, gpoints2, w)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   if (missing(w)) w <- rep(1, n)

   out <- .Fortran("lbtwod",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(b1),as.double(b2),
           as.integer(M1),as.integer(M2), as.double(w), double(M1*M2), PACKAGE="ks")

   xi <- matrix(out[[10]], nrow=M1, ncol=M2)

   return(xi)
}

linbin3D.ksf <- function(X,gpoints1,gpoints2,gpoints3, w)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2],X[,3]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3) 
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   if (missing(w)) w <- rep(1, n)

   out <- .Fortran("lbthrd",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(b1),
           as.double(b2),as.double(b3),as.integer(M1),as.integer(M2),
           as.integer(M3), as.double(w), double(M1*M2*M3),PACKAGE="ks")
   return(array(out[[13]],c(M1,M2,M3)))
}


linbin4D.ksf  <- function(X,gpoints1,gpoints2,gpoints3,gpoints4, w)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2],X[,3],X[,4]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   M3 <- length(gpoints3)
   M4 <- length(gpoints4)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   a3 <- gpoints3[1]
   a4 <- gpoints4[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   b3 <- gpoints3[M3]
   b4 <- gpoints4[M4]
   if (missing(w)) w <- rep(1, n)
   
   out <- .Fortran("lbfoud",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(a4),
           as.double(b1),as.double(b2),as.double(b3),as.double(b4),
           as.integer(M1),as.integer(M2),as.integer(M3),as.integer(M4),
           as.double(w),double(M1*M2*M3*M4),PACKAGE="ks")
   return(array(out[[16]],c(M1,M2,M3,M4)))
}


########################################################################
## Discrete convolution
########################################################################


## Computes the discrete convolution of
## a symmetric or skew-symmetric response 
## vector r and a data vector s.
## If r is symmetric then "skewflag"=1.
## If r is skew-symmetric then "skewflag"=-1.

 
symconv.ks <- function (rr,ss,skewflag = 1) 
{
  L <- length(rr) - 1
  M <- length(ss)
  P <- 2^(ceiling(log(M + L)/log(2)))
  rp <- rep(0,P)
  rp[1:(L+1)] <- rr
  if (L>0) rp[(P-L+1):P] <- skewflag*rr[(L+1):2]
  sp <- rep(0,P) 
  sp[1:M] <- ss
  R <- fft(rp)
  S <- fft(sp)
  t <- fft(R * S, TRUE)
  return((Re(t)/P)[1:M])
}


symconv2D.ks <- function(rr, ss, skewflag=rep(1,2))
{  
  L <- dim(rr)-1
  M <- dim(ss) 
  L1 <- L[1]
  L2 <- L[2]               # find dimensions of r,s
  M1 <- M[1]
  M2 <- M[2]
  P1 <- 2^(ceiling(log(M1+L1)/log(2))) # smallest power of 2 >= M1+L1         
  P2 <- 2^(ceiling(log(M2+L2)/log(2))) # smallest power of 2 >= M2+L2         

  rp <- matrix(0,P1,P2)
  rp[1:(L1+1),1:(L2+1)] <- rr
  if (L1>0)
    rp[(P1-L1+1):P1,1:(L2+1)] <- skewflag[1]*rr[(L1+1):2,]
  if (L2>0)
    rp[1:(L1+1),(P2-L2+1):P2] <- skewflag[2]*rr[,(L2+1):2]
  if (L1 > 0 & L2 > 0)
    rp[(P1-L1+1):P1,(P2-L2+1):P2] <- prod(skewflag)*rr[(L1+1):2,(L2+1):2]   
                                      # wrap around version of rr
  sp <- matrix(0,P1,P2)
  sp[1:M1,1:M2] <- ss                 # zero-padded version of ss

  RR <- fft(rp)        # Obtain FFT's of rr and ss  
  SS <- fft(sp) 
  tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
  return((Re(tt)/(P1*P2))[1:M1,1:M2]) # return normalized truncated tt
}


symconv3D.ks <- function(rr, ss, skewflag=rep(1,3))
{  
   L <- dim(rr) - 1
   M <- dim(ss) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] ; L3 <- L[3]
   M1 <- M[1] ; M2 <- M[2] ; M3 <- M[3]               
   P1 <- P[1] ; P2 <- P[2] ; P3 <- P[3]
   sf <- skewflag

   rp <- array(0,P) 
   rp[1:(L1+1),1:(L2+1),1:(L3+1)] <- rr
   if (L1>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1)]
   if (L2>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1)]
   if (L3>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2]
   if (L1>0 & L2>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1)] <- sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1)]
   if (L2>0 & L3>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3] <- sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2]
   if (L1>0 & L3>0)
     rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3] <- sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2]
   if (L1>0 & L2>0 & L3>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3] <- sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2]

   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3))[1:M1,1:M2,1:M3]) # return normalized truncated tt
}

 
symconv4D.ks <- function(rr, ss, skewflag=rep(1,4) , fftflag=rep(TRUE,2))
{  
   L <- dim(rr) - 1
   M <- dim(ss) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2] ; L3 <- L[3] ; L4 <- L[4]
   M1 <- M[1] ; M2 <- M[2] ; M3 <- M[3] ; M4 <- M[4]               
   P1 <- P[1] ; P2 <- P[2] ; P3 <- P[3] ; P4 <- P[4] 
   sf <- skewflag

   rp <- array(0,P) 
   rp[1:(L1+1),1:(L2+1),1:(L3+1),1:(L4+1)] <- rr

   if (L1>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),1:(L4+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1),1:(L4+1)]
   if (L2>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1),1:(L4+1)]
   if (L3>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2,1:(L4+1)]
   if (L4>0)
     rp[1:(L1+1),1:(L2+1),1:(L3+1),(P4-L4+1):P4] <- sf[4]*rr[1:(L1+1),1:(L2+1),1:(L3+1),(L4+1):2]

   if (L1>0 & L2 >0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <- sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1),1:(L4+1)]
   if (L2>0 & L3>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <- sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2,1:(L4+1)]
   if (L3>0 & L4>0)
     rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,(P4-L4+1):P4] <- sf[3]*sf[4]*rr[1:(L1+1),1:(L2+1),(L3+1):2,(L4+1):2]
   if (L1>0 & L3>0)
     rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <- sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2,1:(L4+1)]
   if (L2>0 & L4>0)
     rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <- sf[2]*sf[4]*rr[1:(L1+1),(L2+1):2,1:(L3+1),(L4+1):2]
   if (L1>0 & L4>0)
     rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),(P4-L4+1):P4] <- sf[1]*sf[4]*rr[(L1+1):2,1:(L2+1),1:(L3+1),(L4+1):2]
   
   if (L1>0 & L2>0 & L3>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <- sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2,1:(L4+1)]
   if (L1>0 & L2>0 & L4>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <- sf[1]*sf[2]*sf[4]*rr[(L1+1):2,(L2+1):2,1:(L3+1),(L4+1):2]
   if (L2>0 & L3>0 & L4>0)
     rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <- sf[2]*sf[3]*sf[4]*rr[1:(L1+1),(L2+1):2,(L3+1):2,(L4+1):2]

   if (L1>0 & L2>0 & L3>0 & L4>0)
     rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <- sf[1]*sf[2]*sf[3]*sf[4]*rr[(L1+1):2,(L2+1):2,(L3+1):2,(L4+1):2]
   
   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3,1:M4] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3*P4))[1:M1,1:M2,1:M3,1:M4]) # return normalized truncated tt
}




