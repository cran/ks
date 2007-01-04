######### R-function:dfltCounts  #########
 
# Obtain default set of grid counts from a 
# multivariate point cloud 'x'.

# Last changed: 18 JUL 2005

dfltCounts.ks <- function(x,gridsize=rep(64,NCOL(x)),h=rep(0,NCOL(x)), supp=1.5, range.x)
{
   x <- as.matrix(x)
   d <- ncol(x)

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

   if (d==1)
      gcounts <- linbin.ks(x,gpoints[[1]])

   if (d==2)
      gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]])

   if (d==3)
      gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]])
   
   if (d==4)
      gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]])
   
   return(list(counts=gcounts,range.x=range.x))
}

######## End of dfltCounts ########

########## R-function: drvkde ##########

# Computes the mth derivative of a binned
# d-variate kernel density estimate based
# on grid counts.

# Last changed: 28 OCT 2005

drvkde.ks <- function(x,drv,bandwidth,gridsize,range.x,binned=FALSE,se=TRUE,estimate.positive=FALSE)
{  
   d <- length(drv)

   if (d==1) x <- as.matrix(x)

   if (binned)
   {
      if (any(dim(x)==1))
         gridsize <- max(dim(x))  

      if (!any(dim(x)==1))
         gridsize <- dim(x)
      gcounts <- x
   }

   if ((!binned)&(missing(gridsize))) gridsize <- rep(64,d)

   # Rename common variables

   h <- bandwidth
   tau <- 4 + max(drv)    

   if (length(h)==1) h <- rep(h,d)

   if (missing(range.x)&(binned==FALSE)) 
   {
      range.x <- list()
      for (id in 1:d)
         range.x[[id]] <- c(min(x[,id])-1.5*h[id],max(x[,id])+1.5*h[id])  
   }

   a <- unlist(lapply(range.x,min))
   b <- unlist(lapply(range.x,max))

   M <- gridsize
   gpoints <- list()

   for (id in 1:d)
      gpoints[[id]] <- seq(a[id],b[id],length=M[id])

   # Bin the data if not already binned

   if (binned==FALSE)
   {
      if (d==1) 
        gcounts <- linbin.ks(x,gpoints[[1]])
      if (d==2) 
        gcounts <- linbin2D.ks(x,gpoints[[1]],gpoints[[2]])
      if (d==3) 
        gcounts <- linbin3D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]])
      if (d==4)
        gcounts <- linbin4D.ks(x,gpoints[[1]],gpoints[[2]],gpoints[[3]],gpoints[[4]])
   }
   else
      gcounts <- x      

   n <- sum(gcounts)
  
   kapmid <- list()
   for (id in (1:d))
   {
     Lid <- min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),M[id]) 
     lvecid <- (0:Lid)
     facid  <- (b[id]-a[id])/(h[id]*(M[id]-1))
     argid <- lvecid*facid
     kapmid[[id]] <- dnorm(argid)/(h[id]^(drv[id]+1))
     hmold0 <- 1
     hmold1 <- argid
     if (drv[id]==0) hmnew <- 1
     if (drv[id]==1) hmnew <- argid
     if (drv[id] >= 2) 
       for (ihm in (2:drv[id])) 
       {
         hmnew <- argid*hmold1 - (ihm-1)*hmold0
         hmold0 <- hmold1   # Compute drv[id] degree Hermite polynomial
         hmold1 <- hmnew    # by recurrence.
       }
     kapmid[[id]] <- hmnew*kapmid[[id]]*(-1)^drv[id]
   }
   
   
   if (d==1)
     kappam <- kapmid[[1]]/n
   
   if (d==2)
     kappam <- outer(kapmid[[1]],kapmid[[2]])/n
     
   if (d==3)
     kappam <- outer(kapmid[[1]],outer(kapmid[[2]],kapmid[[3]]))/n
   

   if (d==4)
     kappam <- outer(kapmid[[1]],
                     outer(kapmid[[2]],outer(kapmid[[3]],kapmid[[4]])))/n
  
   if (!any(c(d==1,d==2,d==3,d==4))) stop("only for d=1,2,3,4")

   if (d==1) 
   { 
      kappam <- as.vector(kappam)
      est <- symconv(kappam,gcounts,skewflag=(-1)^drv)
      
      if (se)
        est.var <- ((symconv((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }

   if (d==2) 
   {     
     est <- symconv2D.ks(kappam,gcounts,skewflag=(-1)^drv)

     if (se)   
       est.var <- ((symconv2D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)     
   }
     

   if (d==3)
   {
     est <- symconv3D.ks(kappam,gcounts,skewflag=(-1)^drv) 
 
     if (se)
       est.var <- ((symconv3D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1)    
   }
     
   if (d==4)
   {
     est <- symconv4D.ks(kappam,gcounts,skewflag=(-1)^drv) 

     if (se)
       est.var <- ((symconv4D.ks((n*kappam)^2,gcounts)/n) - est^2)/(n-1) 
   }
   
   if (estimate.positive)
     est[est<0] <- 0
   
   if (se)
   {
     est.var[est.var<0] <- 0
     return(list(x.grid=gpoints,est=est,se=sqrt(est.var)))
   }
   else if (!se)
     return(list(x.grid=gpoints,est=est))
}

########## End of drvkde #########

########## R-function: linbin ##########

# For application of linear binning to a 
# univariate data set.

# Last changed: 16 JUNE 1995

linbin.ks <- function(X,gpoints,truncate=TRUE)

{
   n <- length(X)
   M <- length(gpoints)  
   trun <- 0
   if (truncate) trun <- 1
   a <- gpoints[1]
   b <- gpoints[M]
   out <- .Fortran("linbin",as.double(X),as.integer(n),
           as.double(a),as.double(b),as.integer(M),
           as.integer(trun),double(M),PACKAGE="ks")
   return(out[[7]])
}

########## End of linbin ##########


######### R-function: linbin2D #########
 
# Creates the grid counts from a bivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbtwod" is called. 

# Last changed: 25 AUG 1995

linbin2D.ks <- function(X,gpoints1,gpoints2)
{
   n <- nrow(X)
   X <- c(X[,1],X[,2]) 
   M1 <- length(gpoints1)
   M2 <- length(gpoints2)
   a1 <- gpoints1[1]
   a2 <- gpoints2[1]
   b1 <- gpoints1[M1]
   b2 <- gpoints2[M2]
   out <- .Fortran("lbtwod",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(b1),as.double(b2),
           as.integer(M1),as.integer(M2),double(M1*M2), PACKAGE="ks")
   return(matrix(out[[9]],M1,M2))
}

########## End of linbin2D ##########

######### R-function: linbin3D #########
 
# Creates the grid counts from a trivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbthrd" is called. 

# Last changed: 27 JUL 2005

linbin3D.ks <- function(X,gpoints1,gpoints2,gpoints3)
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
   out <- .Fortran("lbthrd",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(b1),
           as.double(b2),as.double(b3),as.integer(M1),as.integer(M2),
           as.integer(M3),double(M1*M2*M3),PACKAGE="ks")
   return(array(out[[12]],c(M1,M2,M3)))
}

########## End of linbin3D ##########

######### R-function: linbin4D #########
 
# Creates the grid counts from a quadrivariate data set X 
# over an equally-spaced set of grid points
# contained in "gpoints" using the linear 
# binning strategy. Note that the FORTRAN subroutine
# "lbfoud" is called. 

# Last changed: 31 AUG 2005

linbin4D.ks  <- function(X,gpoints1,gpoints2,gpoints3,gpoints4)
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
   out <- .Fortran("lbfoud",as.double(X),as.integer(n),
           as.double(a1),as.double(a2),as.double(a3),as.double(a4),
           as.double(b1),as.double(b2),as.double(b3),as.double(b4),
           as.integer(M1),as.integer(M2),as.integer(M3),as.integer(M4),
           double(M1*M2*M3*M4),PACKAGE="ks")
   return(array(out[[15]],c(M1,M2,M3,M4)))
}

########## End of linbin4D ##########
#########################################################################
### Binning functions courtesy of M Wand 2005
### Extended by T Duong to 3- and 4-dim 2006
#########################################################################


#
########## R-function: symconv ##########

# Computes the discrete convolution of
# a symmetric or skew-symmetric response 
# vector r and a data vector s.
# If r is symmetric then "skewflag"=1.
# If r is skew-symmetric then "skewflag"=-1.

# Last changed: 03 AUG 2005
 
symconv.ks  <- function(r,s,skewflag=1)

{ 
   L <- length(r)-1
   M <- length(s) 
   P <- 2^(ceiling(log(M+L)/log(2))) # smallest power of 2>=M+L         
   r <- c(r,rep(0,P-2*L-1),skewflag*r[(L+1):2])
                                     # wrap-around version of r
   s <- c(s,rep(0,P-M))              # zero-padded version of s
   R <- fft(r)                       # Obtain FFT's of r and s  
   S <- fft(s)
   t <- fft(R*S,TRUE)               # invert element-wise product of FFT's 
   return((Re(t)/P)[1:M])            # return normalized truncated t
}

########## End of symconv ##########

########## R-function: symconv2D ##########

# Computes the discrete two-dimensional convolution of
# a symmetric response matrix rr and data matrix ss.

# Last changed: 20 MAY 2005
 
symconv2D.ks <- function(rr,ss,skewflag=rep(1,2))

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
   {
      rp[1:(L1+1),(P2-L2+1):P2] <- skewflag[2]*rr[,(L2+1):2]
      rp[(P1-L1+1):P1,(P2-L2+1):P2] <-  prod(skewflag)*rr[(L1+1):2,(L2+1):2]   
   }
                                       # wrap around version of rr
   sp <- matrix(0,P1,P2)
   sp[1:M1,1:M2] <- ss                 # zero-padded version of ss
                                       
   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2))[1:M1,1:M2]) # return normalized truncated tt
}

######## End of symconv2D ########

########### R-function: symconv3D ##########

# Computes the discrete three-dimensional convolution of a symmetric  
# response array rr and a data matrix ss.

# Last changed: 01 JUN 2005
 
symconv3D.ks <- function(rr,ss,skewflag=rep(1,3))

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
   rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1)]
   rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1)]
   rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2]
   rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1)] <-  
                                   sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1)]
   rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3] <-  
                                   sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2]
   rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3] <-  
                                   sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2]
   rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3] <-  
                             sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2]

   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3))[1:M1,1:M2,1:M3]) # return normalized truncated tt
}

########## End of symconv3D ###########

########### R-function: symconv4D ##########

# Computes the discrete four-dimensional convolution of a symmetric  
# response array rr and a data matrix ss.

# Last changed: 01 SEP 2005
 
symconv4D.ks <- function(rr,ss,skewflag=rep(1,4))

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

   rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),1:(L4+1)] <- sf[1]*rr[(L1+1):2,1:(L2+1),1:(L3+1),1:(L4+1)]
   rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <- sf[2]*rr[1:(L1+1),(L2+1):2,1:(L3+1),1:(L4+1)]
   rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <- sf[3]*rr[1:(L1+1),1:(L2+1),(L3+1):2,1:(L4+1)]
   rp[1:(L1+1),1:(L2+1),1:(L3+1),(P4-L4+1):P4] <- sf[4]*rr[1:(L1+1),1:(L2+1),1:(L3+1),(L4+1):2]
     
   rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),1:(L4+1)] <-  
     sf[1]*sf[2]*rr[(L1+1):2,(L2+1):2,1:(L3+1),1:(L4+1)]
   rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <-  
     sf[2]*sf[3]*rr[1:(L1+1),(L2+1):2,(L3+1):2,1:(L4+1)]
   rp[1:(L1+1),1:(L2+1),(P3-L3+1):P3,(P4-L4+1):P4] <-
     sf[3]*sf[4]*rr[1:(L1+1),1:(L2+1),(L3+1):2,(L4+1):2]
   rp[(P1-L1+1):P1,1:(L2+1),(P3-L3+1):P3,1:(L4+1)] <-  
     sf[1]*sf[3]*rr[(L1+1):2,1:(L2+1),(L3+1):2,1:(L4+1)]
   rp[1:(L1+1),(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <-
     sf[2]*sf[4]*rr[1:(L1+1),(L2+1):2,1:(L3+1),(L4+1):2]
   rp[(P1-L1+1):P1,1:(L2+1),1:(L3+1),(P4-L4+1):P4] <-
     sf[1]*sf[4]*rr[(L1+1):2,1:(L2+1),1:(L3+1),(L4+1):2]
   
   rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,1:(L4+1)] <-  
     sf[1]*sf[2]*sf[3]*rr[(L1+1):2,(L2+1):2,(L3+1):2,1:(L4+1)]
   rp[(P1-L1+1):P1,(P2-L2+1):P2,1:(L3+1),(P4-L4+1):P4] <-  
     sf[1]*sf[2]*sf[4]*rr[(L1+1):2,(L2+1):2,1:(L3+1),(L4+1):2]
   rp[1:(L1+1),(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <-  
     sf[2]*sf[3]*sf[4]*rr[1:(L1+1),(L2+1):2,(L3+1):2,(L4+1):2]

   rp[(P1-L1+1):P1,(P2-L2+1):P2,(P3-L3+1):P3,(P4-L4+1):P4] <-  
     sf[1]*sf[2]*sf[3]*sf[4]*rr[(L1+1):2,(L2+1):2,(L3+1):2,(L4+1):2]
   
   sp <- array(0,P)
   sp[1:M1,1:M2,1:M3,1:M4] <- ss            # zero-padded version of ss

   RR <- fft(rp)                       # Obtain FFT's of rr and ss  
   SS <- fft(sp)   
   tt <- fft(RR*SS,TRUE)               # invert element-wise product of FFT's 
   return((Re(tt)/(P1*P2*P3*P4))[1:M1,1:M2,1:M3,1:M4]) # return normalized truncated tt
}

########## End of symconv4D ###########
