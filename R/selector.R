##############################################################################
## Density derivative (psi) functional estimators 
##############################################################################


RKfun <- function(r)
{
  if (r==0)
    val <- 1/(2*sqrt(pi))
  else if (r==1)
    val <- 1/(4*sqrt(pi))
  else if (r==2)
    val <- 3/(8*sqrt(pi)) 
  else if (r==3)
    val <- 15/(16*sqrt(pi))
  else if (r==4)
    val <- 105/(32*sqrt(pi))
  else if (r==5)
    val <- 945/(64*sqrt(pi))
  else if (r==6)
    val <- 10395/(128*sqrt(pi))
  else if (r==7)
    val <- 135135/(256*sqrt(pi))
  else if (r==8)
    val <- 2027025/(512*sqrt(pi))
  
  return(val)
}

#############################################################################
## Estimation of psi_r (no binning, arbitrary d) for G = g^2 I
##
## Code by Jose E. Chacon. Received  04/09/2007
#############################################################################
    
differences <- function(x, upper=TRUE)
{
   n <- nrow(x)
   d <- ncol(x)
   
   difs <- matrix(ncol=d,nrow=n^2)
   for (j in 1:d)
   {    
     xj <- x[,j]
     difxj <- as.vector(xj%*%t(rep(1,n))-rep(1,n)%*%t(xj))
     ##The jth column of difs contains all the differences X_{ij}-X_{kj}
     difs[,j]<-difxj
   }

   if (upper)
   {
     ind.remove <- numeric()
     for (j in 1:(n-1))
       ind.remove <- c(ind.remove, (j*n+1):(j*n+j))

     return(difs[-ind.remove,])
   }
   else
     return(difs)
}

psir.hat <- function(x, r, g, diff=FALSE, upper=diff)
{    
  if(!diff)
  {
    n <- nrow(x)
    d <- ncol(x)
    if (d != length(r))
      stop("The length of r must equal the number of columns of x")
    
    difs <- differences(x, upper=upper)}
  else
  {
    if (upper)
      n <- (-1 + sqrt(1+8*nrow(x)))/2
    else
      n <- sqrt(nrow(x))
    d <- ncol(x)
    difs <- x
  }
  
  drv <- r
  sdrv <- sum(drv)
  arg <- difs/g
  darg <- dmvnorm(arg)/(g^(sdrv+d))
  for (j in 1:d)
  {
    hmold0 <- 1
    hmold1 <- arg[,j]
    hmnew <- 1
    if (drv[j] ==1){hmnew<-hmold1}
    if (drv[j] >= 2) ## Multiply by the corresponding Hermite polynomial, coordinate-wise, using Fact C.1.4 in W&J (1995) and Willink (2005, p.273)
      for (i in (2:drv[j]))
      {
        hmnew <- arg[,j] * hmold1 - (i - 1) * hmold0
        hmold0 <- hmold1
        hmold1 <- hmnew
      }
    darg <- hmnew * darg
  }

  if (diff)
    est <- (sum(darg)*2 - n*darg[1])/n^2
  else
    est <-mean(darg)*(-1)^sdrv
      
  return(est)
}


########################################################################
### Identifying elements of Psi_4 matrix
########################################################################

Psi4.elem <- function(k, kprime, d)
{
  ind <- function(k, d)  
  {
    j <- 1
    dprime <- 1/2*d*(d+1)
    if (k > dprime) stop ("k is larger than d'")
    while (j < d & !(((j-1)*d -1/2*(j-2)*(j-1) < k) & (k <= j*d -1/2*j*(j-1))))
      j <- j+1
    i <- k - (j-1)*d + 1/2*j*(j-1)
    
    return(c(i,j))
  }

  ij <- ind(k, d)
  ei <- elem(ij[1],d)
  ej <- elem(ij[2],d)
  ipjp <- ind(kprime, d)
  eip <- elem(ipjp[1],d)
  ejp <- elem(ipjp[2],d)
  psi4.ind <- ei + eip + ej + ejp
  ##psi4.ind.txt <- paste(psi4.ind, sep="", collapse="")
  coeff <- (2 - (ij[1]==ij[2])) * (2 - (ipjp[1]==ipjp[2]))
  
  return (c(coeff, psi4.ind))
}


Psi4.list <- function(d)
{
  coeff <- vector()
  psifun <- vector()
  dprime <- 1/2*d*(d+1)
  for (k in 1:dprime)
    for (kprime in 1:dprime)
    {
      coeff <- c(coeff, Psi4.elem(k, kprime, d)[1])
      psifun <- rbind(psifun, Psi4.elem(k, kprime, d)[-1]) 
    }

  return(list(coeff=coeff, psi=psifun))  
}


###############################################################################
## Estimate g_AMSE pilot bandwidths for even orders - 2-dim
#
## Parameters
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
#
## Returns
## g_AMSE pilot bandwidths for even orders
###############################################################################

gamse.even.2d <- function(r, n, psi1, psi2)
{
  d <- 2
  num <- -2 * dmvnorm.deriv.2d(x=c(0,0), r=r, Sigma=diag(c(1,1)))
  den <- (psi1 + psi2) * n   
  g.amse <- (num/den)^(1/(2 + d + sum(r)))
  
  return(g.amse)
}

###############################################################################
## Estimate g_AMSE pilot bandwidths for odd orders - 2-dim
#
## Parameters
#
## r - (r1, r2) partial derivative
## n - sample size
## psi1 - psi_(r + (2,0))
## psi2 - psi_(r + (0,2))
## psi00 - psi_(0,0)
## RK - R(K^(r))
## 
## Returns
## g_AMSE pilot bandwidths for odd orders
###############################################################################

gamse.odd.2d <- function(r, n, psi1, psi2, psi00, RK)
{  
    d <- 2
    num <- 2 * psi00 * (2 * sum(r) + d) * RK
    den <- (psi1 + psi2)^2 * n^2
    g.amse <- (num/den)^(1/(2*sum(r) + d + 4))

    return(g.amse)
}

###############################################################################
# Estimate g_SAMSE pilot bandwidth - 2-dim
#
# Parameters
# Sigma.star - scaled variance matrix
# n - sample size
#
# Returns
# g_SAMSE pilot bandwidth
###############################################################################

gsamse.2d <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- 2

  # 4th order derivatives
  if (modr == 4)
  {
    ##RK04 <- 105/(64*pi)
    ##RK13 <- 15/(64*pi)
    ##RK22 <- 9/(64*pi)
    ##RK31 <- 15/(64*pi)
    ##RK40 <- 105/(64*pi)
    ##RK <- c(RK04, RK13, RK22, RK31, RK40)
    
    K40 <- dmvnorm.deriv.2d(x=c(0,0), r=c(4,0), Sigma=diag(c(1,1)))
    K13 <- dmvnorm.deriv.2d(x=c(0,0), r=c(1,3), Sigma=diag(c(1,1)))
    K22 <- dmvnorm.deriv.2d(x=c(0,0), r=c(2,2), Sigma=diag(c(1,1)))
    K31 <- dmvnorm.deriv.2d(x=c(0,0), r=c(3,1), Sigma=diag(c(1,1)))
    K04 <- dmvnorm.deriv.2d(x=c(0,0), r=c(0,4), Sigma=diag(c(1,1)))
    K <- c(K04, K13, K22, K31, K40)
    
    ##psi00 <- psins.2d(r=c(0,0), Sigma=Sigma.star)
    
    if (nstage == 1)
    {
      psi06 <- psins.2d(r=c(0,6), Sigma=Sigma.star)
      psi15 <- psins.2d(r=c(1,5), Sigma=Sigma.star)
      psi24 <- psins.2d(r=c(2,4), Sigma=Sigma.star)
      psi33 <- psins.2d(r=c(3,3), Sigma=Sigma.star)
      psi42 <- psins.2d(r=c(4,2), Sigma=Sigma.star)
      psi51 <- psins.2d(r=c(5,1), Sigma=Sigma.star)
      psi60 <- psins.2d(r=c(6,0), Sigma=Sigma.star)
    }
    else if (nstage == 2)
    {
      psi60 <- psihat[1]
      psi51 <- psihat[2]
      psi42 <- psihat[3]
      psi33 <- psihat[4]
      psi24 <- psihat[5]
      psi15 <- psihat[6]
      psi06 <- psihat[7]
    }
    # required psi functionals    
    psi <- c(psi06 + psi24, psi15 + psi33, psi24 + psi42, psi33 + psi51, 
             psi42 + psi60)
  }
  # 6th order functionals 
  else if (modr == 6)
  {
    ##RK06 <- 10395/(256*pi)
    ##RK15 <- 945/(256*pi)
    ##RK24 <- 315/(256*pi)
    ##RK33 <- 225/(256*pi)
    ##RK42 <- 315/(256*pi)
    ##RK51 <- 945/(256*pi)
    ##RK60 <- 10395/(256*pi)
    ##RK <- c(RK06, RK15, RK24, RK33, RK42, RK51, RK60)
    
    K06 <- dmvnorm.deriv.2d(x=c(0,0), r=c(0,6), Sigma=diag(c(1,1)))
    K15 <- dmvnorm.deriv.2d(x=c(0,0), r=c(1,5), Sigma=diag(c(1,1)))
    K24 <- dmvnorm.deriv.2d(x=c(0,0), r=c(2,4), Sigma=diag(c(1,1)))
    K33 <- dmvnorm.deriv.2d(x=c(0,0), r=c(3,3), Sigma=diag(c(1,1)))
    K42 <- dmvnorm.deriv.2d(x=c(0,0), r=c(4,2), Sigma=diag(c(1,1)))
    K51 <- dmvnorm.deriv.2d(x=c(0,0), r=c(5,1), Sigma=diag(c(1,1)))
    K60 <- dmvnorm.deriv.2d(x=c(0,0), r=c(6,0), Sigma=diag(c(1,1)))
    K <- c(K06, K15, K24, K33, K42, K51, K60)
    
    ##psi00 <- psins.2d(r=c(0,0), Sigma=Sigma.star)
    psi08 <- psins.2d(r=c(0,8), Sigma=Sigma.star)
    psi17 <- psins.2d(r=c(1,7), Sigma=Sigma.star)
    psi26 <- psins.2d(r=c(2,6), Sigma=Sigma.star)
    psi35 <- psins.2d(r=c(3,5), Sigma=Sigma.star)
    psi44 <- psins.2d(r=c(4,4), Sigma=Sigma.star)
    psi53 <- psins.2d(r=c(5,3), Sigma=Sigma.star)
    psi62 <- psins.2d(r=c(6,2), Sigma=Sigma.star)
    psi71 <- psins.2d(r=c(7,1), Sigma=Sigma.star)
    psi80 <- psins.2d(r=c(8,0), Sigma=Sigma.star)

    # required psi functionals
    psi <- c(psi08 + psi26, psi17 + psi35, psi26 + psi44, psi53 + psi35, 
             psi62 + psi44, psi71 + psi53, psi62 + psi80)
  }
 
  ## see thesis for formula
  ##A1 <- sum(RK)*psi00
  A2 <- sum(K^2)
  A3 <- sum(K * psi)  
  A4 <- sum(psi^2)
  B2 <- (2*modr + 2*d)*A2
  B3 <- (modr + d - 2)*A3
  B4 <- A4
  gamma <- (-B3 + sqrt(B3^2 + 4*B2*B4)) / (2*B2)
  
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}



###############################################################################
# Estimate g_SAMSE pilot bandwidth - 3-dim 
#
# Parameters
# Sigma.star - scaled variance matrix
# n - sample size
#
# Returns
# g_SAMSE pilot bandwidth
###############################################################################

gsamse.3d <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- 3
  K <- numeric(); psi <- numeric()

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
          if (sum(c(j1,j2,j3))==d0) derivt <- rbind(derivt, c(j1,j2,j3))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
          if (sum(c(j1,j2,j3))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3))
  d2 <- 8
  derivt8 <- numeric()
  for (j1 in d2:0)
    for (j2 in d2:0)
      for (j3 in d2:0)
          if (sum(c(j1,j2,j3))==d2) derivt8 <- rbind(derivt8, c(j1,j2,j3))

 
                          
  # 4th order g_SAMSE
  if (modr == 4)
  {
    for (i in 1:nrow(derivt))
    {
      r <- derivt[i,]        
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.3d(x=rep(0,d), r=r, Sigma=diag(d)))

        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
            A3psi <- A3psi + psins.3d(r=r+2*elem(j,d), Sigma=Sigma.star)
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          psi <- c(psi, A3psi)     
        }
      }
    }
  }
  #6-th order g_SAMSE
  else if (modr==6)
  {
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        
      
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.3d(x=rep(0,d), r=r, Sigma=diag(d)))
        
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins.3d(r=r+2*elem(j,d), Sigma=Sigma.star)
        
        psi <- c(psi, A3psi)
      }
    }
  }
  
  ##psi00 <- psins.2d(r=rep(0,d), Sigma=Sigma.star)
  ## see thesis for formula
  ## A0 <- sum(RK)
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}


###############################################################################
# Estimate g_SAMSE pilot bandwidth - 4-dim 
#
# Parameters
# Sigma.star - scaled variance matrix
# n - sample size
#
# Returns
# g_SAMSE pilot bandwidth
###############################################################################

gsamse.4d <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- 4
  K <- numeric(); psi <- numeric()

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          if (sum(c(j1,j2,j3,j4))==d0) derivt <- rbind(derivt, c(j1,j2,j3,j4))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          if (sum(c(j1,j2,j3,j4))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3,j4))
  d2 <- 8
  derivt8 <- numeric()
  for (j1 in d2:0)
    for (j2 in d2:0)
      for (j3 in d2:0)
        for (j4 in d2:0)
          if (sum(c(j1,j2,j3,j4))==d2) derivt8 <- rbind(derivt8, c(j1,j2,j3,j4))

 
                          
  ## 4th order g_SAMSE
  if (modr == 4)
  {
    for (i in 1:nrow(derivt))
    {
      r <- derivt[i,]        
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.4d(x=rep(0,d), r=r, Sigma=diag(d)))

        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
            A3psi <- A3psi + psins.4d(r=r+2*elem(j,d), Sigma=Sigma.star)
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          psi <- c(psi, A3psi)     
        }
      }
    }
  }
  ## 6-th order g_SAMSE
  else if (modr==6)
  {
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        
      
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.4d(x=rep(0,d), r=r, Sigma=diag(d)))
        
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins.4d(r=r+2*elem(j,d), Sigma=Sigma.star)
        
        psi <- c(psi, A3psi)
      }
    }
  }
  
  
  ## see thesis for formula
  ## A0 <- sum(RK)
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}


###############################################################################
# Estimate g_SAMSE pilot bandwidth - 5-dim 
#
# Parameters
# Sigma.star - scaled variance matrix
# n - sample size
#
# Returns
# g_SAMSE pilot bandwidth
###############################################################################

gsamse.5d <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- 5
  K <- numeric(); psi <- numeric()

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
              if (sum(c(j1,j2,j3,j4,j5))==d0)
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5))
  d2 <- 8
  derivt8 <- numeric()
  for (j1 in d2:0)
    for (j2 in d2:0)
      for (j3 in d2:0)
        for (j4 in d2:0)
          for (j5 in d2:0)
              if (sum(c(j1,j2,j3,j4,j5))==d2)
                derivt8 <- rbind(derivt8, c(j1,j2,j3,j4,j5))
                        
  ## 4th order g_SAMSE
  if (modr == 4)
  {
    for (i in 1:nrow(derivt))
    {
      r <- derivt[i,]        
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.5d(x=rep(0,d), r=r, Sigma=diag(d)))

        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
            A3psi <- A3psi + psins.5d(r=r+2*elem(j,d), Sigma=Sigma.star)
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          psi <- c(psi, A3psi)
        }
      }
    }
  }
  ## 6-th order g_SAMSE
  else if (modr==6)
  {
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        

      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.5d(x=rep(0,d), r=r, Sigma=diag(d)))
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins.5d(r=r+2*elem(j,d), Sigma=Sigma.star)
        
        psi <- c(psi, A3psi)
      }
    }
  }
  
  ## see thesis for formula
  ## A0 <- sum(RK)
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}



###############################################################################
# Estimate g_SAMSE pilot bandwidth - 6-dim 
#
# Parameters
# Sigma.star - scaled variance matrix
# n - sample size
#
# Returns
# g_SAMSE pilot bandwidth
###############################################################################

gsamse.6d <- function(Sigma.star, n, modr, nstage=1, psihat=NULL)
{
  d <- 6
  K <- numeric(); psi <- numeric()

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
            for (j6 in d0:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d0)
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5,j6))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
            for (j6 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5,j6))
  d2 <- 8
  derivt8 <- numeric()
  for (j1 in d2:0)
    for (j2 in d2:0)
      for (j3 in d2:0)
        for (j4 in d2:0)
          for (j5 in d2:0)
            for (j6 in d2:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d2)
                derivt8 <- rbind(derivt8, c(j1,j2,j3,j4,j5,j6))
                        
  ## 4th order g_SAMSE
  if (modr == 4)
  {
    for (i in 1:nrow(derivt))
    {
      r <- derivt[i,]        
      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.6d(x=rep(0,d), r=r, Sigma=diag(d)))

        A3psi <- 0
        for (j in 1:d)
        {
          if (nstage==1)
            A3psi <- A3psi + psins.6d(r=r+2*elem(j,d), Sigma=Sigma.star)
          else if (nstage==2)
            A3psi <- A3psi + psihat[which.mat(r=r+2*elem(j,d), mat=derivt6)]
          psi <- c(psi, A3psi)
        }
      }
    }
  }
  ## 6-th order g_SAMSE
  else if (modr==6)
  {
    for (i in 1:nrow(derivt6))
    {
      r <- derivt6[i,]        

      if (is.even(r))
      {
        K <- c(K, dmvnorm.deriv.6d(x=rep(0,d), r=r, Sigma=diag(d)))
        A3psi <- 0
        for (j in 1:d)
          A3psi <- A3psi + psins.6d(r=r+2*elem(j,d), Sigma=Sigma.star)
        
        psi <- c(psi, A3psi)
      }
    }
  }
  
  ## see thesis for formula
  ## A0 <- sum(RK)
  A1 <- sum(K^2)
  A2 <- sum(K * psi)  
  A3 <- sum(psi^2)
  B1 <- (2*modr + 2*d)*A1
  B2 <- (modr + d - 2)*A2
  B3 <- A3
  gamma <- (-B2 + sqrt(B2^2 + 4*B1*B3)) / (2*B1)
  
  g.samse <- (gamma * n)^(-1/(modr + d + 2))
  
  return (g.samse)      
}



##############################################################################
# Estimate psi functionals for bivariate data using 1-stage plug-in - 2-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
#
# Returns
# estimated psi functionals
###############################################################################

psifun1.2d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 2
  derivt <- cbind((2*d) - 0:(2*d), 0:(2*d)) 
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  ##RK40 <- 105/(64*pi)
  RK31 <- 15/(64*pi)
  ##RK22 <- 9/(64*pi)
  psi00 <- psins.2d(r=c(0,0), Sigma=S.star) 
  psihat.star <- vector()
  g.star <- vector()

  ## pilots are based on 4th order derivatives
  
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- rep(gsamse.2d(S.star, n, 4), nrow(derivt))

  ## compute 5 different pilots for AMSE
  else if (pilot=="amse")
    for (k in 1:nrow(derivt))
    { 
      r <- derivt[k,]
      psi1 <- psins.2d(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins.2d(r=r + 2*elem(2, 2), Sigma=S.star)

      ## odd order
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK31)
      
      ## even order
      else
        g.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }
  
  if (!binned)
    x.star.diff <- differences(x.star)
   
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    G.star <- g.star[k]^2 * diag(c(1,1))
    if (binned) 
       psihat.star[k] <- dmvnorm.deriv.2d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else
       psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star[k], diff=TRUE)
  }
 
  return(psihat.star)
}




###############################################################################
# Estimate psi functionals for bivariate data using 2-stage plug-in - 2-dim
#
# Parameters
# x - pre-transformed data points
# pilot - "amse" - different AMSE pilot
#       - "samse" - SAMSE pilot
# Returns
# estimated psi functionals
###############################################################################

psifun2.2d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 2
  derivt6 <- cbind((3*d) - 0:(3*d), 0:(3*d)) 
  derivt <- cbind((2*d) - 0:(2*d), 0:(2*d))  
  S.star <- var(x.star)
  n <- nrow(x.star)
  
  ##RK40 <- 105/(64*pi)
  RK31 <- 15/(64*pi)
  ##RK22 <- 9/(64*pi)
  ##RK60 <- 10395/(256*pi)
  RK51 <- 945/(256*pi)
  ##RK42 <- 315/(256*pi)
  RK33 <- 225/(256*pi)
  psi00 <- psins.2d(r=c(0,0), Sigma=S.star) 

  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  g.star <- vector()

  if (!binned)
    x.star.diff <- differences(x.star)
  
  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- rep(gsamse.2d(S.star, n, 6), nrow(derivt6))  
    
  ## compute different pilots for AMSE
  else if (pilot=="amse")
  {       
    for (k in 1:nrow(derivt6))
    {
      r <- derivt6[k,]
      psi1 <- psins.2d(r=r + 2*elem(1, 2), Sigma=S.star)
      psi2 <- psins.2d(r=r + 2*elem(2, 2), Sigma=S.star)
      if (prod(r) == 5)
        g6.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK51)
      else if (prod(r) == 9)
        g6.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK33) 
      else  
        g6.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }
  }
 
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star[k]^2 * diag(c(1,1))
    
    if (binned)
      psihat6.star[k] <- dmvnorm.deriv.2d.sum(x.star, r=r, Sigma=G6.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else
      psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star[k], diff=TRUE)
  }

  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- rep(gsamse.2d(S.star, n, 4, nstage=2,
                              psihat=psihat6.star), nrow(derivt))  
  else if (pilot=="amse")
    for (k in 1:nrow(derivt))
    {
      r <- derivt[k,]
      psi1 <- psihat6.star[7 - (r + 2*elem(1,2))[1]]
      psi2 <- psihat6.star[7 - (r + 2*elem(2,2))[1]]
      
      if (prod(r) == 3)
        g.star[k] <- gamse.odd.2d(r, n, psi1, psi2, psi00, RK31)
      else
        g.star[k] <- gamse.even.2d(r, n, psi1, psi2)
    }

 
 
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    G.star <- g.star[k]^2 * diag(c(1,1))
    if (binned)
      psihat.star[k] <- dmvnorm.deriv.2d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else 
      psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star[k], diff=TRUE)
  }

  return(psihat.star)
}


###############################################################################
# Estimate psi functionals for 3-variate data using 1-stage plug-in - 3-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "samse" = optimal SAMSE pilot bandwidth
# Returns
# estimated psi functionals
###############################################################################

psifun1.3d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 3
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)

  psihat.star <- vector()
  g.star <- vector()

  if (!binned)
    x.star.diff <- differences(x.star)
   
  ## compute 1 pilot for SAMSE
  g.star <- gsamse.3d(S.star, n, 4, nstage=1)
  G.star <- g.star^2 * diag(d)
 
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    if (binned)
      psihat.star[k] <- dmvnorm.deriv.3d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else  
      psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  return(psihat.star)
}



###############################################################################
# Estimate psi functionals for 3-variate data using 2-stage plug-in - 3-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
#
# Returns
# estimated psi functionals
###############################################################################

psifun2.3d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 3

  d0 <- 4
  derivt <- numeric()
    for (j1 in d0:0)
      for (j2 in d0:0)
        for (j3 in d0:0)
            if (sum(c(j1,j2,j3))==d0) derivt <- rbind(derivt, c(j1,j2,j3))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
          if (sum(c(j1,j2,j3))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3))

  S.star <- var(x.star)
  n <- nrow(x.star)

  if (!binned)
    x.star.diff <- differences(x.star)

  psihat6.star <- vector()
  g6.star <- vector()
  psihat.list.star <- vector()
  psihat.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.3d(S.star, n, 6)   

  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
    
    if (binned)
      psihat6.star[k] <- dmvnorm.deriv.3d.sum(x.star, r=r, Sigma=G6.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else
      psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.3d(S.star, n, 4, nstage=2, psihat=psihat6.star) 
  
  ## map psi functionals into the correct order for Psi_4 matrix
  
 
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    G.star <- g.star^2 * diag(d)
    
    if (binned)
     psihat.list.star[k] <- dmvnorm.deriv.3d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else
      psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }  

  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    i <- which.mat(derivt.mat[k,], derivt)
    psihat.star[k] <- psihat.list.star[i]
  }
  
  return(psihat.star)
}



###############################################################################
## Estimate psi functionals for 4-variate data using 1-stage plug-in - 4-dim
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
## Returns
## estimated psi functionals
###############################################################################

psifun1.4d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 4
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)

  if (!binned)
    x.star.diff <- differences(x.star)
  
  psihat.star <- vector()
  g.star <- vector()


  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- gsamse.4d(S.star, n, 4, nstage=1)
  G.star <- g.star^2 * diag(d)
  
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    if (binned)
      psihat.star[k] <- dmvnorm.deriv.4d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else 
      psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)    
  }
  
  return(psihat.star)
}

###############################################################################
# Estimate psi functionals for 4-variate data using 2-stage plug-in - 4-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
#
# Returns
# estimated psi functionals
###############################################################################

psifun2.4d <- function(x.star, pilot="samse", binned, bin.par)
{ 
  d <- 4
  derivt <- numeric()
    for (j1 in d:0)
      for (j2 in d:0)
        for (j3 in d:0)
          for (j4 in d:0)
            if (sum(c(j1,j2,j3,j4))==d) derivt <- rbind(derivt, c(j1,j2,j3,j4))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          if (sum(c(j1,j2,j3,j4))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3,j4))

  S.star <- var(x.star)
  n <- nrow(x.star)

  if (!binned)
    x.star.diff <- differences(x.star)
 
  psihat6.star <- vector()
  g6.star <- vector()
  psihat.list.star <- vector()
  psihat.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
 
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.4d(S.star, n, 6)   
  
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
   
    if (binned)
      psihat.star[k] <- dmvnorm.deriv.4d.sum(x.star, r=r, Sigma=G6.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
    else
      psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.4d(S.star, n, 4, nstage=2, psihat=psihat6.star) 

  
  ## map psi functionals into the correct order for Psi_4 matrix
  
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    G.star <- g.star^2 * diag(d)
  
  if (binned)
    psihat.list.star[k] <- dmvnorm.deriv.3d.sum(x.star, r=r, Sigma=G.star, inc=1, binned=TRUE, bin.par=bin.par)/n^2
  else   
    psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    i <- which.mat(derivt.mat[k,], derivt)
    psihat.star[k] <- psihat.list.star[i]
  }
  
  return(psihat.star)
}

###############################################################################
# Estimate psi functionals for 5-variate data using 1-stage plug-in - 5-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
# Returns
# estimated psi functionals
###############################################################################

psifun1.5d <- function(x.star, pilot="samse")
{ 
  d <- 5
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)

  psihat.star <- vector()
  g.star <- vector()

  x.star.diff <- differences(x.star)
  
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- gsamse.5d(S.star, n, 4, nstage=1)
  ##G.star <- g.star^2 * diag(d)
  
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  return(psihat.star)
}

### psifun1.diag.5d is more efficient than psifun1.5d for diagonal pilot selectors

psifun1.diag.5d <- function(x.star, pilot="samse")
{ 
  d <- 5
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
                           
  psihat.star <- vector()
  psihat.list.star <- vector()
  g.star <- vector()
  

  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- gsamse.5d(S.star, n, 4, nstage=1)
  ##G.star <- g.star^2 * diag(d)

  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    if (is.even(derivt[k,]))
      psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    if (is.even(derivt.mat[k,]))
    {  
      i <- which.mat(derivt.mat[k,], derivt)
      psihat.star[k] <- psihat.list.star[i]
    }
    else
      psihat.star[k] <- 0
  }
  
  return(psihat.star)
}


###############################################################################
# Estimate psi functionals for 5-variate data using 2-stage plug-in - 5-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
#
# Returns
# estimated psi functionals
###############################################################################

psifun2.5d <- function(x.star, pilot="samse")
{ 
  d <- 5
  ##RK <- numeric(); K <- numeric(); psi <- numeric(); 

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
              if (sum(c(j1,j2,j3,j4,j5))==d0)
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5))
 
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  psihat.list.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.5d(S.star, n, 6)   

  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }    
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.5d(S.star, n, 4, nstage=2, psihat=psihat6.star) 
 
  
  ## map psi functionals into the correct order for Psi_4 matrix     
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    i <- which.mat(derivt.mat[k,], derivt)
    psihat.star[k] <- psihat.list.star[i]
  }
  
  return(psihat.star)
}



### psifun2.diag.5d is more efficient than psifun2.5d for diagonal pilot selectors

psifun2.diag.5d <- function(x.star, pilot="samse")
{ 
  d <- 5
  ##RK <- numeric(); K <- numeric(); psi <- numeric(); 

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
              if ((sum(c(j1,j2,j3,j4,j5))==d0) && is.even(c(j1,j2,j3,j4,j5)))
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5))
 
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  psihat.list.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.5d(S.star, n, 6)   

  psihat6.star <- rep(0, nrow(derivt6))
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    if (is.even(r))
      psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }    
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.5d(S.star, n, 4, nstage=2, psihat=psihat6.star) 
 
  
  ## map psi functionals into the correct order for Psi_4 matrix     
  for (k in 1:nrow(derivt))
  { 
    r <- derivt[k,] 
    if (is.even(r))
      psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    if (is.even(derivt.mat[k,]))
    {  
      i <- which.mat(derivt.mat[k,], derivt)
      psihat.star[k] <- psihat.list.star[i]
    }
    else
      psihat.star[k] <- 0
  }
  
  return(psihat.star)
}

###############################################################################
# Estimate psi functionals for 6-variate data using 1-stage plug-in - 6-dim
#
# Parameters
# x.star - pre-transformed data points
# pilot - "amse" = different AMSE pilot bandwidths
#       - "samse" = optimal SAMSE pilot bandwidth
# Returns
# estimated psi functionals
###############################################################################

psifun1.6d <- function(x.star, pilot="samse")
{ 
  d <- 6
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat.star <- vector()
  g.star <- vector()

  
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- gsamse.6d(S.star, n, 4, nstage=1)
  ##G.star <- g.star^2 * diag(d)
  
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,]
    psihat.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  return(psihat.star)
}

### psifun1.diag.6d is more efficient than psifun1.6d for diagonal pilot selectors

psifun1.diag.6d <- function(x.star, pilot="samse")
{ 
  d <- 6
  derivt <- Psi4.list(d)$psi
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat.star <- vector()
  psihat.list.star <- vector()

 
  ## compute 1 pilot for SAMSE
  if (pilot=="samse")
    g.star <- gsamse.6d(S.star, n, 4, nstage=1)
  ##G.star <- g.star^2 * diag(d)
  
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    if (is.even(derivt[k,]))
      psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    if (is.even(derivt.mat[k,]))
    {  
      i <- which.mat(derivt.mat[k,], derivt)
      psihat.star[k] <- psihat.list.star[i]
    }
    else
      psihat.star[k] <- 0
    
  }
  
  return(psihat.star)
}


###############################################################################
## Estimate psi functionals for 6-variate data using 2-stage plug-in - 6-dim
##
## Parameters
## x.star - pre-transformed data points
## pilot - "amse" = different AMSE pilot bandwidths
##       - "samse" = optimal SAMSE pilot bandwidth
##
## Returns
## estimated psi functionals
###############################################################################

psifun2.6d <- function(x.star, pilot="samse")
{ 
  d <- 6
  ##RK <- numeric(); K <- numeric(); psi <- numeric(); 

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
            for (j6 in d0:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d0)
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5,j6))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
            for (j6 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5,j6))
 
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  psihat.list.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  ## compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.6d(S.star, n, 6)   

  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }    
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.6d(S.star, n, 4, nstage=2, psihat=psihat6.star) 
 
  
  ## map psi functionals into the correct order for Psi_4 matrix     
  for (k in 1:nrow(derivt))
  {
    r <- derivt[k,] 
    psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    i <- which.mat(derivt.mat[k,], derivt)
    psihat.star[k] <- psihat.list.star[i]
  }
  
  return(psihat.star)
}



### psifun2.diag.6d is more efficient than psifun2.6d for diagonal selectors

psifun2.diag.6d <- function(x.star, pilot="samse")
{ 
  d <- 6
  ##RK <- numeric(); K <- numeric(); psi <- numeric(); 

  d0 <- 4
  derivt <- numeric()
  for (j1 in d0:0)
    for (j2 in d0:0)
      for (j3 in d0:0)
        for (j4 in d0:0)
          for (j5 in d0:0)
            for (j6 in d0:0)
              if ((sum(c(j1,j2,j3,j4,j5,j6))==d0) && is.even(c(j1,j2,j3,j4,j5,j6)))
                derivt <- rbind(derivt, c(j1,j2,j3,j4,j5,j6))
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
            for (j6 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5,j6))
 
  S.star <- var(x.star)
  n <- nrow(x.star)
  x.star.diff <- differences(x.star)
  
  psihat6.star <- vector()
  g6.star <- vector()
  psihat.star <- vector()
  psihat.list.star <- vector()
  g.star <- vector()

  ## pilots are based on 6th order derivatives
  
  # compute 1 pilot for SAMSE    
  if (pilot=="samse")
    g6.star <- gsamse.6d(S.star, n, 6)   

  psihat6.star <- rep(0, nrow(derivt6))
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    if (is.even(r))
      psihat6.star[k] <- psir.hat(x=x.star.diff, r=r, g=g6.star, diff=TRUE)
  }    
  
  ## pilots are based on 4th order derivatives using 6th order psi functionals
  ## computed above 'psihat6.star'
    
  if (pilot=="samse")
    g.star <- gsamse.6d(S.star, n, 4, nstage=2, psihat=psihat6.star) 
 
  
  ## map psi functionals into the correct order for Psi_4 matrix     
  for (k in 1:nrow(derivt))
  { 
    r <- derivt[k,] 
    if (is.even(r))
      psihat.list.star[k] <- psir.hat(x=x.star.diff, r=r, g=g.star, diff=TRUE)
    
  }
  
  derivt.mat <- Psi4.list(d)$psi
  for (k in 1:nrow(derivt.mat))
  {
    if (is.even(derivt.mat[k,]))
    {  
      i <- which.mat(derivt.mat[k,], derivt)
      psihat.star[k] <- psihat.list.star[i]
    }
    else
      psihat.star[k] <- 0
  }
  
  return(psihat.star)
}

###############################################################################
## Estimate psi functionals for 6-variate data using 1-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
###############################################################################


psifun1.unconstr.nd <- function(x, Sd4, Sd6, rel.tol=10^-10)
{
  n <- nrow(x)
  d <- ncol(x)
  S <- var(x)
  
  nlim <- 1e4
  upper <- TRUE 

  ## stage 1 of plug-in
  G4 <-(2^(d/2+3)/((d+4)*n))^(2/(d+8))*S
  vecPsi4 <- vecPsir(x=x, Sdr=Sd4, Gr=G4, r=4, upper=upper, nlim=nlim)
  
  return (vecPsi4)
}


###############################################################################
## Estimate psi functionals for 6-variate data using 2-stage plug-in 
## with unconstrained pilot
##
## Parameters
## x - data points
## Sd4, Sd6 - symmetrizer matrices of order 4 and 6
##
## Returns
## estimated psi functionals
###############################################################################


psifun2.unconstr.nd <- function(x, Sd4, Sd6, rel.tol=10^-10)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  Hstart <- (4/(d+2))^(2/(d+4))*n^(-2/(d+4))*S
  Hstart <- matrix.sqrt(Hstart)
  nlim <- 1e4
 
  ## matrix of pairwise differences
  upper <- TRUE
  difs <- differences(x, upper=upper)
 
  ## constants for normal reference
  D4phi0 <- D4L0(d=d, Sd4=Sd4)
  Id1 <- diag(d)
  vId <- vec(Id1)

  ## stage 1 of plug-in
  G6 <- (2^(d/2+5)/((d+6)*n))^(2/(d+8))*S
  G612 <- matrix.sqrt(G6)
 
  vecPsi6 <- vecPsir(x=x, Sdr=Sd6, Gr=G6, r=6, upper=upper, nlim=nlim)
   
  Id4 <- diag(d^4)
  Id2 <- diag(d^2)
  Kdd2 <- K.mat(m=d,n=d^2)
  Psi6 <- (Id1%x%Kdd2%x%Id2)%*%vecPsi6
  
  ## asymptotic squared bias for r = 4
  AB2r4<-function(vechG){
    r <- 4
    G <- invvech(vechG)%*%invvech(vechG)
    G12 <- matrix.sqrt(G)
    Ginv12 <- chol2inv(chol(G12))
    AB <- n^(-1)*det(Ginv12)*(K.pow(A=Ginv12,pow=r)%*%D4phi0)+
      (1/2)*(t(vec(G))%x%Id4)%*%Psi6
    return (sum(AB^2))
  }
  
  res <- optim(vech(Hstart),AB2r4, control=list(reltol=rel.tol))
  V4 <- res$value
  G4 <- res$par
  G4 <- invvech(G4)%*%invvech(G4)
 
  ## stage 2 of plug-in

  vecPsi4 <- vecPsir(x=x, Sdr=Sd4, Gr=G4, r=4, upper=upper, nlim=nlim)
  
  return (vecPsi4)
}

### vec Psi_r functionals estimation with unconstr bw matrices

vecPsir <- function(x, Gr, Sdr, r, upper, nlim=1e4)
{
  d <- ncol(x)
  n <- nrow(x)
  S <- var(x)

  ## matrix of differences - upper triangular form is available
  difs <- differences(x, upper=upper)
  
  Id1 <- diag(d)
  vId <- vec(Id1)
  Gr12 <- matrix.sqrt(Gr)
  Grinv12 <- chol2inv(chol(Gr12))

  args <- t(Grinv12%*%t(difs))

  if (r==6)
  {  
    vecPsi6 <- rep(0,d^6)

    K6 <- function(args)
    {
      return(mat.K.pow(args,6)-15*mat.Kprod(mat.K.pow(args,4),rep(1,nrow(args))%x%t(vId)) + 45*mat.Kprod(mat.K.pow(args,2),rep(1,nrow(args))%x%t(K.pow(vId,2)))-15*rep(1,nrow(args))%x%t(K.pow(vId,3)))
    }

    ## split into blocks because of memory limitations
    if (nrow(args) > nlim)
    {
      num.loops <- nrow(args) %/% nlim
      for (j in 1:num.loops)
      {
        args.temp <- args[((j-1)* nlim+1):(j*nlim),]
        hmnew.temp <- K6(args.temp)
        vecPsi6  <- vecPsi6 + colSums(dmvnorm(args.temp)*hmnew.temp)
        ##cat(j, " ")
      }
    
      if (nrow(args) %% nlim >0)
      {
        args.temp <- args[(num.loops*nlim+1):nrow(args),]
        hmnew.temp <- K6(args.temp)
        vecPsi6  <- vecPsi6 + colSums(dmvnorm(args.temp)*hmnew.temp)
      }
    }
    else
    {
      hmnew <- K6(args)
      vecPsi6 <- colSums(dmvnorm(args)*hmnew)
    }
    ##cat("\n")

    if (upper)
    {
      ## adjust for upper triangular form of differences matrix
      args0 <- t(rep(0,d))
      hmnew0 <- K6(args0)
      vecPsi6 <- 2*vecPsi6 - as.vector(n*dmvnorm(args0)*hmnew0)
      vecPsi6 <- Sdr%*%vecPsi6/n^2
    }
    else
      vecPsi6 <- Sdr%*%vecPsi6/nrow(args)

    vecPsi6 <- (-1)^6*det(Grinv12)*K.pow(Grinv12,6)%*%vecPsi6

    return (vecPsi6)
  }

  
  if (r==4)
  {
    vecPsi4 <- rep(0,d^4)
    K4 <- function(args)
    {
      return(mat.K.pow(args,4)-6*mat.Kprod(mat.K.pow(args,2),rep(1,nrow(args))%x%t(vId))+3*rep(1,nrow(args))%x%t(K.pow(vId,2)))
    }
    
    if (nrow(args)>nlim)
    {
      num.loops <- nrow(args) %/% nlim
      for (j in 1:num.loops)
      {
        args.temp <- args[((j-1)* nlim+1):(j*nlim),]
        hmnew.temp <- K4(args.temp)
        vecPsi4  <- vecPsi4 + colSums(dmvnorm(args.temp)*hmnew.temp)
        ##cat(j, " ")
      }
    
      if (nrow(args) %% nlim >0)
      {
        args.temp <- args[(num.loops*nlim+1):nrow(args),]
        hmnew.temp <- K4(args.temp)
        vecPsi4  <- vecPsi4 + colSums(dmvnorm(args.temp)*hmnew.temp)
      }   
    }
    else
    {
      hmnew <- K4(args)
      vecPsi4 <- colSums(dmvnorm(args)*hmnew)
    }
    ##cat("\n")
    
    if (upper)
    {
      args0 <- t(rep(0,d))
      hmnew0 <- K4(args0)
      vecPsi4 <- 2*vecPsi4 - as.vector(n*dmvnorm(args0)*hmnew0)
      vecPsi4 <- Sdr%*%vecPsi4/n^2
    }
    else
      vecPsi4 <- Sdr%*%vecPsi4/nrow(args)

    vecPsi4<-(-1)^4*det(Grinv12)*K.pow(Grinv12,4)%*%vecPsi4

    return(vecPsi4)
  }
}


###############################################################################
## Creates Psi_4 matrix of 4th order psi functionals used in AMISE - 2 to 6 dim
##
## Parameters
## x - data points
## nstage - number of plug-in stages (1 or 2)
## pilot - "amse" - different AMSE pilot
##       - "samse" - SAMSE pilot
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## matrix of psi functionals
############################################################################

psimat.2d <- function(x.star, nstage=1, pilot="samse", binned, bin.par)
{
  if (nstage==1)
    psi.fun <- psifun1.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  else if (nstage==2)
    psi.fun <- psifun2.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  
  psi40 <- psi.fun[1] 
  psi31 <- psi.fun[2] 
  psi22 <- psi.fun[3] 
  psi13 <- psi.fun[4] 
  psi04 <- psi.fun[5] 

  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13, psi04)
  
  return(matrix(coeff * psi.fun, nc=3, nr=3))
}

psimat.3d <- function(x.star, nstage=1, pilot="samse", binned, bin.par)
{
  d <- 3
  
  if (nstage==1)
    psi.fun <- psifun1.3d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  else if (nstage==2)
    psi.fun <- psifun2.3d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  
  coeff <- Psi4.list(3)$coeff
  
  return(matrix(coeff * psi.fun, nc=d*(d+1)/2, nr=d*(d+1)/2))
}

psimat.4d <- function(x.star, nstage=1, pilot="samse", binned, bin.par)
{
  d <- 4 
 
  if (nstage==1)
    psi.fun <- psifun1.4d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  else if (nstage==2)
    psi.fun <- psifun2.4d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
  
  coeff <- Psi4.list(d)$coeff
 
  return(matrix(coeff * psi.fun, nc=d*(d+1)/2, nr=d*(d+1)/2))
}


psimat.5d <- function(x.star, nstage=1, pilot="samse")
{
  d <- 5

  if (nstage==1)
    psi.fun <- psifun1.5d(x.star, pilot=pilot)
  else if (nstage==2)
    psi.fun <- psifun2.5d(x.star, pilot=pilot)
  
  coeff <- Psi4.list(d)$coeff

  return(matrix(coeff * psi.fun,  nc=d*(d+1)/2, nr=d*(d+1)/2))
}

### psimat.diag.5d is more efficient than psimat.5d for diagonal selectors

psimat.diag.5d <- function(x.star, nstage=1, pilot="samse")
{
  d <- 5
 
  if (nstage==1)
    psi.fun <- psifun1.diag.5d(x.star, pilot=pilot)
  else if (nstage==2)
    psi.fun <- psifun2.diag.5d(x.star, pilot=pilot)

  coeff <- Psi4.list(d)$coeff

  return(matrix(coeff * psi.fun, nc=d*(d+1)/2, nr=d*(d+1)/2))
}

psimat.6d <- function(x.star, nstage=1, pilot="samse")
{
  d <- 6
  
  if (nstage==1)
    psi.fun <- psifun1.6d(x.star, pilot=pilot)
  else if (nstage==2)
    psi.fun <- psifun2.6d(x.star, pilot=pilot)
  
  coeff <- Psi4.list(d)$coeff

  return(matrix(coeff * psi.fun,  nc=d*(d+1)/2, nr=d*(d+1)/2))
}

### psimat.diag.6d is more efficient than psimat.6d for diagonal selectors

psimat.diag.6d <- function(x.star, nstage=1, pilot="samse")
{
  d <- 6
 
  if (nstage==1)
    psi.fun <- psifun1.diag.6d(x.star, pilot=pilot)
  else if (nstage==2)
    psi.fun <- psifun2.diag.6d(x.star, pilot=pilot)

  coeff <- Psi4.list(d)$coeff

  return(matrix(coeff * psi.fun, nc=d*(d+1)/2, nr=d*(d+1)/2))
}

###  unconstrained pilot selectors

psimat.unconstr.nd <- function(x, nstage=1, Sd4, Sd6)
{
  if (nstage==1)
    psi.fun <- psifun1.unconstr.nd(x=x, Sd4=Sd4, Sd6=Sd6)
  else if (nstage==2)
    psi.fun <- psifun2.unconstr.nd(x=x, Sd4=Sd4, Sd6=Sd6)

  return(invvec(psi.fun))
}
  
###############################################################################
# Plug-in bandwidth selectors
###############################################################################

    
###############################################################################
## Computes plug-in full bandwidth matrix - 2 to 6 dim
##
## Parameters
## x - data points
## Hstart - initial value for minimisation
## nstage - number of plug-in stages (1 or 2)
## pilot - "amse" - different AMSE pilot
##       - "samse" - SAMSE pilot
##       - "unconstr" - unconstrained pilot
## pre - "scale" - pre-scaled data
##     - "sphere"- pre-sphered data 
##
## Returns
## Plug-in full bandwidth matrix
###############################################################################

hpi <- function(x, nstage=2, binned=TRUE, bgridsize)
{
  ## 1-d selector is taken from KernSmooth's dpik
  
  if (missing(bgridsize)) bgridsize <- 401
  return(dpik(x=x, level=nstage, gridsize=bgridsize))
}

Hpi <- function(x, nstage=2, pilot="samse", pre="sphere", Hstart, binned=FALSE, bgridsize, amise=FALSE)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)

  if(!is.matrix(x)) x <- as.matrix(x)

  if (substr(pre,1,2)=="sc")
  {
    x.star <- pre.scale(x)
    S12 <- diag(sqrt(diag(var(x))))
    Sinv12 <- chol2inv(chol(S12))
  }
  else if (substr(pre,1,2)=="sp")
  {
    x.star <- pre.sphere(x)
    S12 <- matrix.sqrt(var(x))
    Sinv12 <- chol2inv(chol(S12))
  }
  
  if (substr(pilot,1,1)=="a")
    pilot <- "amse"
  else if (substr(pilot,1,1)=="s")
    pilot <- "samse"
  else if (substr(pilot,1,1)=="u")
    pilot <- "unconstr"
           
  
  if (pilot=="amse" & d>2)
    stop("SAMSE pilot selectors are better for higher dimensions")

  if (pilot=="unconstr" & d>=6)
    stop("Uconstrained pilots not implemented yet for 6-dim data")
  
  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  if (d > 4) binned <- FALSE
  
  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    ## linear binning
    bin.par <- dfltCounts.ks(x.star, bgridsize, sqrt(diag(H.max)))
  }

  ## psi4.mat is on pre-transformed data scale
  if (pilot!="unconstr")
  {  
    if (d==2)
      psi4.mat <- psimat.2d(x.star, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par)
    else if (d==3)
      psi4.mat <- psimat.3d(x.star, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par)
    else if (d==4)
      psi4.mat <- psimat.4d(x.star, nstage=nstage, pilot=pilot, binned=binned, bin.par=bin.par)
    else if (d==5)
      psi4.mat <- psimat.5d(x.star, nstage=nstage, pilot=pilot)
    else if (d==6)
      psi4.mat <- psimat.6d(x.star, nstage=nstage, pilot=pilot)
  }
  else
  {
    ##require(Matrix)
    ##data(Sd4)
    ##data(Sd6)
    ##Sd4 <- Sd4[[d]]
    ##Sd6 <- Sd6[[d]]
    ##require(combinat)
    Sd2 <- Sdr(d=d, r=2)
    Sd4 <- Sdr(d=d, r=4)
    Sd6 <- Sdr(d=d, r=6)
    psi4.mat <- psimat.unconstr.nd(x=x, nstage=nstage, Sd4=Sd4, Sd6=Sd6)
  }

  if (pilot=="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x)
    
    Hstart <- matrix.sqrt(Hstart)

    ## PI is estimate of AMISE
    pi.temp.unconstr <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vec(H)) %*% psi4.mat %*% vec(H)
      return(drop(pi.temp))
    }

    ## psi4.mat always a zero eigen-values since it has repeated rows 
    result <- optim(vech(Hstart), pi.temp.unconstr, method="BFGS")
    H <- invvech(result$par) %*% invvech(result$par)
  }
  else if (pilot!="unconstr")
  {
    ## use normal reference bandwidth as initial condition 
    if (missing(Hstart)) 
      Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else
      Hstart <- Sinv12 %*% Hstart %*% Sinv12
    
    Hstart <- matrix.sqrt(Hstart)

    ## PI is estimate of AMISE
    pi.temp <- function(vechH)
    { 
      H <- invvech(vechH) %*% invvech(vechH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vech(H)) %*% psi4.mat %*% vech(H)
      return(drop(pi.temp)) 
    }
    
    ## check that Psi_4 is positive definite (needed for AMSE pilots) 
    if (prod(eigen(psi4.mat)$val > 0) == 1)
    {
      result <- optim(vech(Hstart), pi.temp, method="BFGS")
      H <- invvech(result$par) %*% invvech(result$par)
    }
    else
    { 
      cat("Psi matrix not positive definite\n")
      H <- matrix(NA, nc=d, nr=d) 
    }    
    ## back-transform
    H <- S12 %*% H %*% S12
  }
  
  if (!amise)
    return(H)
  else
    return(list(H = H, PI=result$value))
}     



###############################################################################
# Computes plug-in diagonal bandwidth matrix for 2 to 6-dim
#
# Parameters
# x - data points
# nstage - number of plug-in stages (1 or 2)
# pre - "scale" - pre-scaled data
#     - "sphere"- pre-sphered data 
#
# Returns
# Plug-in diagonal bandwidth matrix
###############################################################################

Hpi.diag <- function(x, nstage=2, pilot="amse", pre="scale", Hstart, binned=FALSE,
                     bgridsize)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  
  if (substr(pre,1,2)=="sc")
    x.star <- pre.scale(x)
  else if (substr(pre,1,2)=="sp")
    x.star <- pre.sphere(x)

  if (substr(pre,1,2)=="sp")
    stop("Using pre-sphering won't give diagonal bandwidth matrix\n")

  if (substr(pilot,1,1)=="a")
    pilot <- "amse"
  else if (substr(pilot,1,1)=="s")
    pilot <- "samse"
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  s1 <- sd(x[,1])
  s2 <- sd(x[,2])

  if (substr(pre,1,2)=="sc") S12 <- diag(sqrt(diag(var(x))))
  else if (substr(pre,1,2)=="sp") S12 <- matrix.sqrt(var(x))
 
  Sinv12 <- chol2inv(chol(S12))
    
  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  if (d > 4) binned <- FALSE
  
  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
    ## linear binning
    bin.par <- dfltCounts.ks(x.star, bgridsize, sqrt(diag(H.max)))
  }
  
  if (d==2)
  {
    if (nstage == 1)
      psi.fun <- psifun1.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
    else if (nstage == 2)
      psi.fun <- psifun2.2d(x.star, pilot=pilot, binned=binned, bin.par=bin.par)
    
    psi40 <- psi.fun[1]
    psi22 <- psi.fun[3]
    psi04 <- psi.fun[5]
    
    ## diagonal bandwidth matrix for 2-dim has exact formula 
    h1 <- (psi04^(3/4)*RK/(psi40^(3/4)*(sqrt(psi40 * psi04)+psi22)*n))^(1/6)
    h2 <- (psi40/psi04)^(1/4) * h1

    return(diag(c(s1^2*h1^2, s2^2*h2^2)))
  }
  else 
  { 
    if (pilot=="amse")
      stop("SAMSE pilot selectors are better for higher dimensions")
 
    ## use normal reference bandwidth as initial condition

    if (missing(Hstart)) 
       Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
    else    
       Hstart <- Sinv12 %*% Hstart %*% Sinv12

    Hstart <- matrix.sqrt(Hstart)
 
    if (d==3)
      psi4.mat <- psimat.3d(x.star, nstage=nstage, pilot=pilot,binned=binned, bin.par=bin.par)    
    else if (d==4)
      psi4.mat <- psimat.4d(x.star, nstage=nstage, pilot=pilot,binned=binned, bin.par=bin.par)
    else if (d==5)
      psi4.mat <- psimat.diag.5d(x.star, nstage=nstage, pilot=pilot)
    else if (d==6)
      psi4.mat <- psimat.diag.6d(x.star, nstage=nstage, pilot=pilot)

   
    ## PI is estimate of AMISE
    pi.temp <- function(diagH)
    { 
      H <- diag(diagH) %*% diag(diagH)
      pi.temp <- 1/(det(H)^(1/2)*n)*RK + 1/4* t(vech(H)) %*% psi4.mat %*% vech(H)
    return(drop(pi.temp)) 
    }
    
    result <- optim(diag(Hstart), pi.temp, method="BFGS")
    H <- diag(result$par) %*% diag(result$par)
  
  ## back-transform
  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))
  H <- S12 %*% H %*% S12
  
  return(H)
  }
}



###############################################################################
# Cross-validation bandwidth selectors
###############################################################################

###############################################################################
# Computes the least squares cross validation LSCV function for 2 to 6 dim
# 
# Parameters
# x - data values
# H - bandwidth matrix
#
# Returns
# LSCV(H)
###############################################################################

lscv.1d.binned <- function(xbin.par, h)
{
  n <- sum(xbin.par$counts)
  lscv1 <- n^2*bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], drv=0, bandwidth=sqrt(2)*h, binned=TRUE)
  lscv2 <- n^2*(bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], drv=0, bandwidth=h, binned=TRUE) - dnorm(0, mean=0, sd=h)/n)
 
  lscv <- lscv1/n^2 - 2/(n*(n-1))*lscv2

  return(lscv)
}



lscv.mat <- function(x, H, binned=FALSE, bin.par)
{
  n <- nrow(x)
  d <- ncol(x)

  ##h <- sqrt(diag(H))

  if (d==2)
  {
    lscv1 <- dmvnorm.2d.sum(x, 2*H, inc=1, binned=binned, bin.par=bin.par)
    lscv2 <- dmvnorm.2d.sum(x, H, inc=0, binned=binned, bin.par=bin.par)
  }
  else if (d==3)
  {
    lscv1 <- dmvnorm.3d.sum(x, 2*H, inc=1, binned=binned, bin.par=bin.par)
    lscv2 <- dmvnorm.3d.sum(x, H, inc=0, binned=binned, bin.par=bin.par)
  }
  else if (d==4)
  {
    lscv1 <- dmvnorm.4d.sum(x, 2*H, inc=1, binned=binned, bin.par=bin.par)
    lscv2 <- dmvnorm.4d.sum(x, H, inc=0, binned=binned, bin.par=bin.par)
  }
  else if (d==5)
  {
    lscv1 <- dmvnorm.5d.sum(x, 2*H, inc=1)
    lscv2 <- dmvnorm.5d.sum(x, H, inc=0)
  }
  else if (d==6)
  {
    lscv1 <- dmvnorm.6d.sum(x, 2*H, inc=1)
    lscv2 <- dmvnorm.6d.sum(x, H, inc=0)
  }
  
  return(lscv1/n^2 - 2/(n*(n-1))*lscv2)     
}

   
###############################################################################
# Finds the bandwidth matrix that minimises LSCV for 2 to 6 dim
# 
# Parameters
# x - data values
# Hstart - initial bandwidth matrix
#
# Returns
# H_LSCV
###############################################################################

hlscv <- function(x, binned=TRUE, bgridsize)
{
  if (any(duplicated(x)))
    stop("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- length(x)
  d <- 1
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))
  if (missing(bgridsize)) bgridsize <- 401
  
  xbin.par <- dfltCounts.ks(x, gridsize=bgridsize)

  lscv.1d.temp <- function(h)
  {
    return(lscv.1d.binned(x=xbin.par, h=h))
  }

  opt <- optimise(f=lscv.1d.temp, interval=c(0.2*hnorm, 5*hnorm, tol=.Machine$double.eps))$minimum
  return(opt)
    
}
  
Hlscv <- function(x, Hstart)
{
  if (any(duplicated(x)))
    stop("Data contain duplicated values: LSCV is not well-behaved in this case")
  n <- nrow(x)
  d <- ncol(x)
  ##RK <- (4*pi)^(-d/2)
  
  ## use normal reference selector as initial condn
  if (missing(Hstart)) 
    Hstart <- matrix.sqrt((4/ (n*(d + 2)))^(2/(d + 4)) * var(x))

  lscv.mat.temp <- function(vechH)
  {
    ##  ensures that H is positive definite
    H <- invvech(vechH) %*% invvech(vechH)
    return(lscv.mat(x=x, H=H, binned=FALSE))
  }
  result <- optim(vech(Hstart), lscv.mat.temp, method="Nelder-Mead")
                                        #control=list(abstol=n^(-10*d)))    
  
  return(invvech(result$par) %*% invvech(result$par))
}

###############################################################################
# Finds the diagonal bandwidth matrix that minimises LSCV for 2 to 6 dim
# 
# Parameters
# x - data values
# Hstart - initial bandwidth matrix
#
# Returns
# H_LSCV,diag
###############################################################################

Hlscv.diag <- function(x, Hstart, binned=FALSE, bgridsize)
{
  n <- nrow(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  if (missing(Hstart)) 
    Hstart <- matrix.sqrt((4/ (n*(d + 2)))^(2/(d + 4)) * var(x))

  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  if (d > 4) binned <- FALSE

  if (binned)
  {
    H.max <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x)
    ## linear binning
    bin.par <- dfltCounts.ks(x, bgridsize, sqrt(diag(H.max)))
  }
  
  lscv.mat.temp <- function(diagH)
  {
    H <- diag(diagH^2)
    ## ensures that H is positive definite

    return(lscv.mat(x, H, binned=binned, bin.par=bin.par))
  }
  result <- optim(diag(Hstart), lscv.mat.temp, method="Nelder-Mead")
                  #control=list(abstol=n^(-10*d)))   
  return(diag(result$par^2))
}

###############################################################################
# Computes the biased cross validation BCV function for 2-dim
# 
# Parameters
# x - data values
# H1, H2 - bandwidth matrices
#
# Returns
# BCV(H)
###############################################################################

bcv.mat <- function(x, H1, H2)
{
  n <- nrow(x)
  d <- 2
  
  psi40 <- dmvnorm.deriv.2d.sum(x, Sigma=H2, r=c(4,0), inc=0)
  psi31 <- dmvnorm.deriv.2d.sum(x, Sigma=H2, r=c(3,1), inc=0)
  psi22 <- dmvnorm.deriv.2d.sum(x, Sigma=H2, r=c(2,2), inc=0)
  psi13 <- dmvnorm.deriv.2d.sum(x, Sigma=H2, r=c(1,3), inc=0)
  psi04 <- dmvnorm.deriv.2d.sum(x, Sigma=H2, r=c(0,4), inc=0)
    
  coeff <- c(1, 2, 1, 2, 4, 2, 1, 2, 1)
  psi.fun <- c(psi40, psi31, psi22, psi31, psi22, psi13, psi22, psi13,psi04)/
    (n*(n-1))
  psi4.mat <- matrix(coeff * psi.fun, nc=3, nr=3)
  
  RK <- (4*pi)^(-d/2) 
  bcv <- drop(n^(-1)*det(H1)^(-1/2)*RK + 1/4*t(vech(H1)) %*% psi4.mat
              %*% vech(H1))
  
  return(list(bcv=bcv, psimat=psi4.mat))
}



###############################################################################
# Find the bandwidth matrix that minimises the BCV for 2-dim
# 
# Parameters
# x - data values
# whichbcv - 1 = BCV1
#          - 2 = BCV2 
# Hstart - initial bandwidth matrix
#
# Returns
# H_BCV
###############################################################################

Hbcv <- function(x, whichbcv=1, Hstart)
{
  n <- nrow(x)
  d <- ncol(x)
  D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  RK <- (4*pi)^(-d/2)

  # use normal reference b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- Hmax
  ##lo.bound <- -Hmax
  
  if (missing(Hstart))
    Hstart <- matrix.sqrt(0.9*Hmax)

  bcv1.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    # ensures that H is positive definite

    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(vechH)
  {
    H <- invvech(vechH) %*% invvech(vechH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }

  # derivatives of BCV1 function - see thesis  
  bcv1.mat.deriv <- function(vechH)
  {
    H <-  invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    
    psi4.mat <- bcv.mat(x, H, H)$psimat
    psi22 <- psi4.mat[1,3] 
    psi00 <- dmvnorm.2d.sum(x, Sigma=H, inc=0)/(n*(n-1))
    psi22.deriv.xxt <- dmvnorm.deriv.2d.xxt.sum(x, r=c(2,2), Sigma=H)/(n*(n-1))
    psi22.deriv <- t(D2)%*% vec((Hinv %*% psi22.deriv.xxt %*% Hinv +
                                 2* psi00 *Hinv %*% Hinv - psi22*Hinv)/2) 
    
    const <- matrix(c(0,0,1, 0,4,0, 1,0,0), nc=3, byrow=TRUE)
    psi4.mat.deriv<- const %x% psi22.deriv
    
    deriv1 <- -1/2*n^{-1}*RK*t(D2) %*% vec(chol2inv(chol(H)))
    deriv2 <- 1/2 * psi4.mat %*% vech(H) + 1/4 *
      (t(psi4.mat.deriv) %*% (vech(H) %x% diag(c(1,1,1)))) %*% vech(H)
    
    return(deriv1 + deriv2)      
  }
  
  # derivatives of BCV2 function - see thesis   
  bcv2.mat.deriv <- function(vechH)
  {
    H <-  invvech(vechH) %*% invvech(vechH)
    Hinv <- chol2inv(chol(H))
    
    psi4.mat <- bcv.mat(x, H, 2*H)$psimat
    psi22 <- psi4.mat[1,3] 
    psi00 <- dmvnorm.2d.sum(x, Sigma=2*H, inc=0)/(n*(n-1))
    psi22.deriv.xxt <- dmvnorm.deriv.2d.xxt.sum(x,r=c(2,2),Sigma=2*H)/(n*(n-1))
    psi22.deriv <- t(D2)%*% vec((Hinv %*% psi22.deriv.xxt %*% Hinv +
                                 2* psi00 *Hinv %*% Hinv - psi22*Hinv)/2) 
    
    const <- matrix(c(0,0,1, 0,4,0, 1,0,0), nc=3, byrow=TRUE)
    psi4.mat.deriv<- const %x% psi22.deriv
    
    deriv1 <- -1/2*n^{-1}*RK*t(D2) %*% vec(chol2inv(chol(H)))
    deriv2 <- 1/2 * psi4.mat %*% vech(H) + 1/4 *
      (t(psi4.mat.deriv) %*% (vech(H) %x% diag(c(1,1,1)))) %*% vech(H)
    
    return(deriv1 + 2*deriv2)    
  }
  
  if (whichbcv==1)
    result <- optim(vech(Hstart), bcv1.mat.temp, gr=bcv1.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)))
  else if (whichbcv==2)
    result <- optim(vech(Hstart), bcv2.mat.temp, gr=bcv2.mat.deriv,
                    method="L-BFGS-B", upper=vech(matrix.sqrt(up.bound)),
                    lower=-vech(matrix.sqrt(up.bound)))
  
  return(invvech(result$par) %*% invvech(result$par))
}

###############################################################################
# Find the diagonal bandwidth matrix that minimises the BCV for 2-dim
# 
# Parameters
# x - data values
# whichbcv - 1 = BCV1
#          - 2 = BCV2
# Hstart - initial bandwidth matrix
#
# Returns
# H_BCV, diag
###############################################################################

Hbcv.diag <- function(x, whichbcv=1, Hstart)
{
  n <- nrow(x)
  d <- ncol(x)
  ##D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  RK <- (4*pi)^(-d/2)
  
  ## use maximally smoothed b/w matrix for bounds
  k <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*n*gamma((d+8)/2)*(d+2)))^(2/(d+4))
  Hmax <- k * abs(var(x))
  up.bound <- diag(Hmax)
  ##lo.bound <- rep(0,d)
  
  if (missing(Hstart))
    Hstart <- 0.9*matrix.sqrt(Hmax)

  bcv1.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    ## ensures that H is positive definite

    return(bcv.mat(x, H, H)$bcv)
  }
    
  bcv2.mat.temp <- function(diagH)
  {
    H <- diag(diagH) %*% diag(diagH)
    return(bcv.mat(x, H, 2*H)$bcv)
  }
  
  if (whichbcv == 1)
    result <- optim(diag(Hstart), bcv1.mat.temp, method="L-BFGS-B",
                    upper=sqrt(up.bound))
  else if (whichbcv == 2)
    result <- optim(diag(Hstart), bcv2.mat.temp, method="L-BFGS-B",
                    upper=sqrt(up.bound))

  return(diag(result$par) %*% diag(result$par))
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
  

###############################################################################
# Estimate g_AMSE pilot bandwidth for SCV for 2 to 6 dim
#
# Parameters
# Sigma.star - scaled/ sphered variance matrix
# Hamise - (estimate) of H_AMISE 
# n - sample size
#
# Returns
# g_AMSE pilot bandwidth
###############################################################################

gamse.scv.2d <- function(x.star, Sigma.star, Hamise, n, binned=FALSE, bgridsize)
{
  d <- 2
  ##derivt6 <- cbind((3*d) - 0:(3*d), 0:(3*d))
  g6.star <- gsamse.2d(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(c(1,1))

  ## required psi functionals

  if (binned)
    bin.par <- dfltCounts.ks(x.star, bgridsize, sqrt(diag(G6.star)))
        
  psi60 <- dmvnorm.deriv.2d.sum(x.star, r=c(6,0), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi51 <- dmvnorm.deriv.2d.sum(x.star, r=c(5,1), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi42 <- dmvnorm.deriv.2d.sum(x.star, r=c(4,2), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi33 <- dmvnorm.deriv.2d.sum(x.star, r=c(3,3), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi24 <- dmvnorm.deriv.2d.sum(x.star, r=c(2,4), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi15 <- dmvnorm.deriv.2d.sum(x.star, r=c(1,5), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  psi06 <- dmvnorm.deriv.2d.sum(x.star, r=c(0,6), Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)/n^2
  
  Theta6 <- invvech(c(psi60 + 2*psi42 + psi24, psi51 + 2*psi33 + psi15,
                     psi42 + 2*psi24 + psi06))
  eye2 <- diag(c(1,1))
  D2 <- rbind(c(1,0,0), c(0,1,0), c(0,1,0), c(0,0,1))
  trHamise <- Hamise[1,1] + Hamise[2,2]

  # required constants - see thesis
  Cmu1 <- 1/2*t(D2) %*% vec(Theta6 %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D2)%*% vec(Hamise)
                               + trHamise * t(D2) %*% vec(eye2))
  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) +
    sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}


gamse.scv.3d <- function(x.star, Sigma.star, Hamise, n, binned=FALSE, bgridsize)
{
  d <- 3
  psi.fun6 <- vector()
  g6.star <- gsamse.3d(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(d)
  
  ## required psi functionals
  if (binned)
    bin.par <- dfltCounts.ks(x.star, bgridsize, sqrt(diag(G6.star)))
  
  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
          if (sum(c(j1,j2,j3))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3))

  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
    psihat6 <- dmvnorm.deriv.3d.sum(x.star, r=r, Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)
    psi.fun6[k] <- psihat6/(n^2)
  }    
  
  Theta6.mat <- matrix(0, nc=d, nr=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
      {
        temp.sum <- temp.sum + psi.fun6[which.mat(temp[k,], derivt6)]
      }
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye3 <- diag(d)
  D4 <- dupl(d)$d
  trHamise <- Hamise[1,1] + Hamise[2,2] + Hamise[3,3] 

  ## required constants - see thesis
  Cmu1 <- 1/2*t(D4) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D4)%*% vec(Hamise)
                               + trHamise * t(D4) %*% vec(eye3))

  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) +
    sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}



gamse.scv.4d <- function(x.star, Sigma.star, Hamise, n, binned=FALSE, bgridsize)
{
  d <- 4
  psi.fun6 <- vector()
  g6.star <- gsamse.4d(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(d)

  if (binned)
    bin.par <- dfltCounts.ks(x.star, bgridsize, sqrt(diag(G6.star)))
  
  ## required psi functionals

  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          if (sum(c(j1,j2,j3,j4))==d1) derivt6 <- rbind(derivt6, c(j1,j2,j3,j4))

  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
    psihat6 <- dmvnorm.deriv.4d.sum(x.star, r=r, Sigma=G6.star, inc=1, binned=binned, bin.par=bin.par)
    psi.fun6[k] <- psihat6/(n^2)
  }    
  
  Theta6.mat <- matrix(0, nc=d, nr=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
      {
        temp.sum <- temp.sum + psi.fun6[which.mat(temp[k,], derivt6)]
      }
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye4 <- diag(d)
  D4 <- dupl(d)$d
  trHamise <- Hamise[1,1] + Hamise[2,2] + Hamise[3,3] + Hamise[4,4]

  # required constants - see thesis
  Cmu1 <- 1/2*t(D4) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D4)%*% vec(Hamise)
                               + trHamise * t(D4) %*% vec(eye4))
  
  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) +
    sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}

gamse.scv.5d <- function(x.star, Sigma.star, Hamise, n)
{
  d <- 5
  g6.star <- gsamse.5d(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(d)
  psi.fun6 <- vector()
  
  # required psi functionals

  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5))
 
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
    psihat6 <- dmvnorm.deriv.5d.sum(x.star, r=r, Sigma=G6.star, inc=1)
    psi.fun6[k] <- psihat6/(n^2)
  }    
     
  Theta6.mat <- matrix(0, nc=d, nr=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
      {
        temp.sum <- temp.sum + psi.fun6[which.mat(temp[k,], derivt6)]
      }
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye5 <- diag(d)
  D6 <- dupl(d)$d
  trHamise <- Hamise[1,1] + Hamise[2,2] + Hamise[3,3] + Hamise[4,4]+ Hamise[5,5]
  
  # required constants - see thesis
  Cmu1 <- 1/2*t(D6) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D6)%*% vec(Hamise)
                               + trHamise * t(D6) %*% vec(eye5))
 
  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) +
    sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}

gamse.scv.6d <- function(x.star, Sigma.star, Hamise, n)
{
  d <- 6
  g6.star <- gsamse.6d(Sigma.star, n, 6) 
  G6.star <- g6.star^2 * diag(d)
  psi.fun6 <- vector()
  
  # required psi functionals

  d1 <- 6
  derivt6 <- numeric()
  for (j1 in d1:0)
    for (j2 in d1:0)
      for (j3 in d1:0)
        for (j4 in d1:0)
          for (j5 in d1:0)
            for (j6 in d1:0)
              if (sum(c(j1,j2,j3,j4,j5,j6))==d1)
                derivt6 <- rbind(derivt6, c(j1,j2,j3,j4,j5,j6))
 
  for (k in 1:nrow(derivt6))
  {
    r <- derivt6[k,]
    G6.star <- g6.star^2 * diag(d)
    psihat6 <- dmvnorm.deriv.6d.sum(x.star, r=r, Sigma=G6.star, inc=1)
    psi.fun6[k] <- psihat6/(n^2)
  }    
     
  Theta6.mat <- matrix(0, nc=d, nr=d)
  Theta6.mat.ind <- Theta6.elem(d)
  for (i in 1:d)
    for (j in 1:d)
    {
      temp <- Theta6.mat.ind[[i]][[j]]
      temp.sum <- 0
      for (k in 1:nrow(temp))
      {
        temp.sum <- temp.sum + psi.fun6[which.mat(temp[k,], derivt6)]
      }
      Theta6.mat[i,j] <- temp.sum 
    }
    
  eye6 <- diag(d)
  D6 <- dupl(d)$d
  trHamise <- Hamise[1,1] + Hamise[2,2] + Hamise[3,3] + Hamise[4,4]+ Hamise[5,5] + Hamise[6,6]
  
  # required constants - see thesis
  Cmu1 <- 1/2*t(D6) %*% vec(Theta6.mat %*% Hamise)
  Cmu2 <- 1/8*(4*pi)^(-d/2) * (2*t(D6)%*% vec(Hamise)
                               + trHamise * t(D6) %*% vec(eye6))

  num <- 2 * (d+4) * sum(Cmu2*Cmu2)
  den <- -(d+2) * sum(Cmu1*Cmu2) +
    sqrt((d+2)^2 * sum(Cmu1*Cmu2)^2 + 8*(d+4)*sum(Cmu1*Cmu1) * sum(Cmu2*Cmu2))
  gamse <- (num / (den*n))^(1/(d+6)) 

  return(gamse)
}


###############################################################################
# Computes the smoothed cross validation function for 2 to 6 dim
# 
# Parameters
# x - data values
# H - bandwidth matrix
# G - pilot bandwidth matrix
#
# Returns
# SCV(H)
###############################################################################


scv.1d <- function(x, h, g, binned=TRUE)
{
  if (binned)
  {  
    n <- sum(x$counts)
    scv1 <- dnorm.1d.sum(bin.par=x, sigma=sqrt(2*h^2+2*g^2), binned=binned)
    scv2 <- dnorm.1d.sum(bin.par=x, sigma=sqrt(h^2+2*g^2), binned=binned)
    scv3 <- dnorm.1d.sum(bin.par=x, sigma=sqrt(2*g^2), binned=binned)
    
    #if (scv1 <0) scv1 <- 1e10
    #if (scv2 <0) scv2 <- 1e10
    #if (scv3 <0) scv3 <- 1e10
  }
  else
  {
    n <- length(x)
    scv1 <- dnorm.1d.sum(x=x, sigma=sqrt(2*h^2+2*g^2), binned=FALSE)
    scv2 <- dnorm.1d.sum(x=x, sigma=sqrt(h^2+2*g^2), binned=FALSE)
    scv3 <- dnorm.1d.sum(x=x, sigma=sqrt(2*g^2), binned=FALSE)
  }
  scv <- n^(-2)*(scv1 - 2*scv2 + scv3) + (n*h)^(-1)*(4*pi)^(-1/2) ##+ 

  return(scv)
}

scv.mat <- function(x, H, G)
{
  n <- nrow(x)
  d <- ncol(x)

  if (d==2)
  {
    scv1 <- dmvnorm.2d.sum(x, Sigma=2*H + 2*G, inc=1)
    scv2 <- dmvnorm.2d.sum(x, Sigma=H + 2*G, inc=1)
    scv3 <- dmvnorm.2d.sum(x, Sigma=2*G, inc=1)
  }
  else if (d==3)
  {

    scv1 <- dmvnorm.3d.sum(x, Sigma=2*H + 2*G, inc=1)
    scv2 <- dmvnorm.3d.sum(x, Sigma=H + 2*G, inc=1)
    scv3 <- dmvnorm.3d.sum(x, Sigma=2*G, inc=1)
  }
  else if (d==4)
  {
    scv1 <- dmvnorm.4d.sum(x, Sigma=2*H + 2*G, inc=1)
    scv2 <- dmvnorm.4d.sum(x, Sigma=H + 2*G, inc=1)
    scv3 <- dmvnorm.4d.sum(x, Sigma=2*G, inc=1)
  } 
  else if (d==5)
  {
    scv1 <- dmvnorm.5d.sum(x, Sigma=2*H + 2*G, inc=1)
    scv2 <- dmvnorm.5d.sum(x, Sigma=H + 2*G, inc=1)
    scv3 <- dmvnorm.5d.sum(x, Sigma=2*G, inc=1)
  } 
  else if (d==6)
  {
    scv1 <- dmvnorm.6d.sum(x, Sigma=2*H + 2*G, inc=1)
    scv2 <- dmvnorm.6d.sum(x, Sigma=H + 2*G, inc=1)
    scv3 <- dmvnorm.6d.sum(x, Sigma=2*G, inc=1)
  }
  scvmat <- n^(-1)*det(H)^(-1/2)*(4*pi)^(-d/2)+ n^(-2)*(scv1 - 2*scv2 + scv3)
    
  return (scvmat)
}



###############################################################################
# Find the bandwidth matrix that minimises the SCV for 2 to 6 dim
# 
# Parameters
# x - data values
# pre - "scale" - pre-scaled data
#     - "sphere"- pre-sphered data
# Hstart - initial bandwidth matrix
#
# Returns
# H_SCV
###############################################################################

hscv <- function(x, nstage=2, binned=TRUE, bgridsize, hmin, hmax)
{
  sigma <- sd(x)
  n <- length(x)
  d <- 1
  hnorm <- sqrt((4/(n*(d + 2)))^(2/(d + 4)) * var(x))
  if (missing(bgridsize)) bgridsize <- 401
  if (missing(hmin)) hmin <- 0.1*hnorm
  if (missing(hmax)) hmax <- 2*hnorm
  
  xbin.par <- dfltCounts.ks(x, gridsize=bgridsize)

  if (nstage==1)
  {
    psihat6 <- psins.1d(r=6, sigma=sigma)
    psihat10 <- psins.1d(r=10, sigma=sigma)
  }
  else if (nstage==2)
  {
    psihat8 <- psins.1d(r=8, sigma=sigma)
    psihat12 <- psins.1d(r=12, sigma=sigma)
    g1 <- (2/(7*n))^(1/9)*2^(1/2)*sigma
    g2 <- (2/(11*n))^(1/13)*2^(1/2)*sigma
    
    psihat6 <- bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], binned=TRUE, drv=6, bandwidth=g1)
    psihat10 <- bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], binned=TRUE, drv=10, bandwidth=g2) 
  }

  g3 <- (-6/((2*pi)^(1/2)*psihat6*n))^(1/7) 
  g4 <- (-210/((2*pi)^(1/2)*psihat10*n))^(1/11)
  psihat4 <- bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], binned=TRUE, drv=4, bandwidth=g3)
  psihat8 <- bkfe(x=xbin.par$counts, range.x=xbin.par$range.x[[1]], binned=TRUE, drv=8, bandwidth=g4)   
  C <- (441/(64*pi))^(1/18) * (4*pi)^(-1/5) * psihat4^(-2/5) * psihat8^(-1/9)
  if (n > 1e5) h0 <- hpi(x)
  
  scv.1d.temp <- function(h)
  {
    ##print(c(h, C*n^(-23/45)*h^(-2)))
    return(scv.1d(x=xbin.par, h=h, g=C*n^(-23/45)*h^(-2), binned=TRUE))
  }

  #plot(0, xlim=c(hmin, hmax), type="n", ylim=c(0,5*scv.1d.temp(hpi(x))))
  #for (i in seq(hmin, hmax, length=100))
  #  points(i, scv.1d.temp(i))
  #browser()
  opt <- optimise(f=scv.1d.temp, interval=c(hmin, hmax))$minimum
  return(opt)
}


Hscv <- function(x, pre="sphere", Hstart, binned=FALSE, bgridsize)
{
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)

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

  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  if (d > 4) binned <- FALSE
 
  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))

  S12inv <- chol2inv(chol(S12))
  Hamise <- S12inv %*% Hpi(x=x,nstage=1,pilot="samse", pre="sphere", binned=binned, bgridsize=bgridsize) %*% S12inv

  if (any(is.na(Hamise)))
  {
    warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
    Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
  }

  if (d==2)
    gamse <- gamse.scv.2d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==3)
    gamse <- gamse.scv.3d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==4)
    gamse <- gamse.scv.4d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==5)
    gamse <- gamse.scv.5d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n)
  else if (d==6)
    gamse <- gamse.scv.6d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n)
  
  G.amse <- gamse^2 * diag(d)
  
  ## use normal reference bandwidth as initial condition
  if (missing(Hstart)) 
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
  else    
    Hstart <- S12inv %*% Hstart %*% S12inv

  Hstart <- matrix.sqrt(Hstart)

  scv.mat.temp <- function(vechH)
  {
    ## ensures that H is positive definite
    H <- invvech(vechH) %*% invvech(vechH)
    return(scv.mat(x.star, H, G.amse))
  }
  
  ## back-transform
  result <- optim(vech(Hstart), scv.mat.temp, method= "Nelder-Mead")
                                        #control=list(abstol=n^(-10*d)))
 
  H <- invvech(result$par) %*% invvech(result$par)
  H <- S12 %*% H %*% S12
  
  
  return(H)
}


Hscv.diag <- function(x, pre="scale", Hstart, binned=FALSE, bgridsize)
{
  if(!is.matrix(x)) x <- as.matrix(x)
  d <- ncol(x)
  RK <- (4*pi)^(-d/2)
  
  ## pre-transform data
  
  if (substr(pre,1,2)=="sc")
    pre <- "scale"
  else if (substr(pre,1,2)=="sp")
    pre <- "sphere"

  if (pre=="sphere")
     stop("Using pre-sphering doesn't give a diagonal bandwidth matrix\n")
  
  if (pre=="sphere")
    x.star <- pre.sphere(x)
  else if (pre=="scale")
    x.star <- pre.scale(x)
  
  S.star <- var(x.star)
  n <- nrow(x.star)

  if (missing(bgridsize) & binned)
    if (d==2)
      bgridsize <- rep(151,d)
    else if (d==3)
      bgridsize <- rep(51, d)
    else if (d==4)
      bgridsize <- rep(21, d)

  if (d > 4) binned <- FALSE
   
  if (pre=="scale") S12 <- diag(sqrt(diag(var(x))))
  else if (pre=="sphere") S12 <- matrix.sqrt(var(x))

  S12inv <- chol2inv(chol(S12))
  Hamise <- S12inv %*% Hpi(x=x,nstage=1,pilot="samse", pre="sphere", binned=binned, bgridsize=bgridsize) %*% S12inv

  if (any(is.na(Hamise)))
  {
    warning("Pilot bandwidth matrix is NA - replaced with maximally smoothed")
    Hamise <- (((d+8)^((d+6)/2)*pi^(d/2)*RK)/(16*(d+2)*n*gamma(d/2+4)))^(2/(d+4))* var(x.star)
  }

  if (d==2)
    gamse <- gamse.scv.2d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==3)
    gamse <- gamse.scv.3d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==4)
    gamse <- gamse.scv.4d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n, binned=binned, bgridsize=bgridsize)
  else if (d==5)
    gamse <- gamse.scv.5d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n)
  else if (d==6)
    gamse <- gamse.scv.6d(x.star=x.star, Sigma.star=S.star, H=Hamise, n=n)
  
  G.amse <- gamse^2 * diag(d)
  
  ## use normal reference bandwidth as initial condition
  if (missing(Hstart)) 
    Hstart <- (4/(n*(d + 2)))^(2/(d + 4)) * var(x.star)
  else    
    Hstart <- S12inv %*% Hstart %*% S12inv

  Hstart <- matrix.sqrt(Hstart)

  scv.mat.temp <- function(diagH)
  {
    ## ensures that H is positive definite
    H <- diag(diagH) %*% diag(diagH)
    return(scv.mat(x.star, H, G.amse))
  }
  
  ## back-transform
  result <- optim(diag(Hstart), scv.mat.temp, method= "Nelder-Mead")
                                        #control=list(abstol=n^(-10*d)))
  H <- diag(result$par) %*% diag(result$par)
  H <- S12 %*% H %*% S12
  
  return(H)
}
