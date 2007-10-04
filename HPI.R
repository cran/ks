library(mvtnorm)
library(ks)

##### Odd factorial and matrix square root

OF<-function(m){factorial(m)/(2^(m/2)*factorial(m/2))} 

matrix.sqrt <- function(A) {
    d<-ncol(A)
    if(d==1){A12<-sqrt(A)}
    else{
    xe <- eigen(A)
    xe1 <- xe$values
    if(all(xe1 >= 0)) {
    xev1 <- diag(sqrt(xe1))
    }
    xval1 <- cbind(xe$vectors)
    xval1i <- solve(xval1)
    A12 <- xval1 %*% xev1 %*% xval1i}
    return(A12)
    }

    
###### Estimation of psi_r (no binning, arbitrary d)
    
psir.hat <- function(x, r, g){
    n<-nrow(x)
    d<-ncol(x)
    if(d!=length(r)){stop("The length of r must equal the number of columns of x")}
    difs<-matrix(ncol=d,nrow=n^2)
    for(j in 1:d){    
        xj<-x[,j]
        difxj<-as.vector(xj%*%t(rep(1,n))-rep(1,n)%*%t(xj)) #The jth column of difs contains all the differences X_{ij}-X_{kj}
        difs[,j]<-difxj}

    drv<-r
    sdrv<-sum(drv)
    arg<-difs/g
    darg<-dmvnorm(arg)/(g^(sdrv+d))
    for(j in 1:d){
        hmold0 <- 1
        hmold1 <- arg[,j]
        hmnew <- 1
        if (drv[j] ==1){hmnew<-hmold1}
        if (drv[j] >= 2) ## Multiply by the corresponding Hermite polynomial, coordinate-wise, using Fact C.1.4 in W&J (1995) and Willink (2005, p.273)
            for (i in (2:drv[j])) {
                hmnew <- arg[,j] * hmold1 - (i - 1) * hmold0
                hmold0 <- hmold1
                hmold1 <- hmnew
                }
    darg <- hmnew * darg
    }
    est<-mean(darg)*(-1)^sdrv
    return(est)}

##### Psir for a multivariate normal density
    
psinorm<-function(r){### General d, Sigma=Id
    d<-length(r)
    sr<-sum(r)
    all.even<-sum(r-2*floor(r/2))
    if(all.even==0){res<-(-1)^(sr/2)*prod(OF(r))*(4*pi)^(-d/2)*2^(-sr/2)}
    else{res<-0}
    return(res)
    }

psinorm.2d<-function(m,Sigma){### d=2, general Sigma, following the formulas in Wand and Jones (1994)
    s1<-sqrt(Sigma[1,1])
    s2<-sqrt(Sigma[2,2])
    rho<-Sigma[1,2]/(s1*s2)
    sm<-sum(m)
    if(2*floor(sm/2)!=sm){psi<-0}
    else{
        if(2*floor(m[1]/2)==m[1]){
            r<-m[1]/2
            s<-m[2]/2
            j<-0:min(r,s)
            lambda<-(factorial(2*r)*factorial(2*s)/2^(r+s))*sum((2*rho)^(2*j)/(factorial(r-j)*factorial(s-j)*factorial(2*j)))
            }
        else{
            r<-(m[1]-1)/2
            s<-(m[2]-1)/2
            j<-0:min(r,s)
            lambda<-(rho*factorial(2*r+1)*factorial(2*s+1)/2^(r+s))*
                sum((2*rho)^(2*j)/(factorial(r-j)*factorial(s-j)*factorial(2*j+1)))
            }
        psi<-(-1)^(sm/2)*lambda/(2^((sm+4)/2)*pi*s1^(m[1]+1)*s2^(m[2]+1)*(1-rho^2)^((sm+1)/2))
        }
    return(psi)
    }

 
###### 2-stage estimators of psi_0 and psi_r (2-dimensional case), with AMSE pilots

psi0.2d<-function(samp){
    d<-ncol(samp)
    n<-nrow(samp)
    S<-var(samp)
    psis<-numeric(d)
    for(i in 1:d){
        z1<-rep(0,length=d)
        z1[i]<-2
        psis.aux<-numeric(d)
            for(j in 1:d){
                z2<-rep(0,length=d)
                z2[j]<-2
                psis.aux[j]<-psinorm.2d(m=z1+z2,Sigma=S)
                }
        spsis.aux<-sum(psis.aux)
        sz1<-2
        Kr0<-(-1)^sz1*psinorm.2d(m=z1,Sigma=diag(2)/2)
        az1<-(-2*Kr0/(spsis.aux*n))^(1/(d+sz1+2))
        psis[i]<-psir.hat(samp,r=z1,g=az1)
        }
    spsis<-sum(psis)
    zeros<-rep(0,length=d)
    sr<-sum(zeros)
    Kr0<-(-1)^sr*psinorm.2d(m=zeros,Sigma=diag(2)/2)
    a<-(-2*Kr0/(spsis*n))^(1/(d+sr+2))
    est<-psir.hat(samp,r=zeros,g=a)
    return(est)
    }
    
psir.2d<-function(samp,r,psi0){
    d<-ncol(samp)
    n<-nrow(samp)
    S<-var(samp)
    if(missing(psi0)){psi0<-psi0.2d(samp)}
    psis<-numeric(d)
    for(i in 1:d){
        z1<-rep(0,length=d)
        z1[i]<-2
        rz1<-r+z1
        psis.aux<-numeric(d)
            for(j in 1:d){
                z2<-rep(0,length=d)
                z2[j]<-2
                psis.aux[j]<-psinorm.2d(m=rz1+z2,Sigma=S)
                }
        spsis.aux<-sum(psis.aux)
        sz1<-sum(rz1)
        Kr0<-(-1)^sz1*psinorm.2d(m=rz1,Sigma=diag(2)/2)
        RKr<-(-1)^sz1*psinorm.2d(m=2*rz1,Sigma=diag(2))
        all.even<-sum(rz1-2*floor(rz1/2))
        if(all.even==0){
            az1<-(-2*Kr0/(spsis.aux*n))^(1/(d+sz1+2))}
        else{
            az1<-((2*psi0*(2*sz1+d)*RKr)/(spsis.aux^2*n^2))^(1/(2*sz1+d+4))}
        psis[i]<-psir.hat(samp,r=rz1,g=az1)
        }
    spsis<-sum(psis)
    sr<-sum(r)
    Kr0<-(-1)^sr*psinorm.2d(m=r,Sigma=diag(2)/2)
    RKr<-(-1)^sr*psinorm.2d(m=2*r,Sigma=diag(2))
    all.even<-sum(r-2*floor(r/2))
    if(all.even==0){
        az2<-(-2*Kr0/(spsis*n))^(1/(d+sr+2))}
    else{
        az2<-((2*psi0*(2*sr+d)*RKr)/(spsis^2*n^2))^(1/(2*sr+d+4))}
    est<-psir.hat(samp,r=r,g=az2)
    return(est)
    }
    

###### Estimation of the matrix Psi4 with AMSE pilots

Psi4<-function(samp){
    psi0<-psi0.2d(samp)
    psi4s<-numeric(5)
    for(i in 0:4){
        psi4s[i+1]<-psir.2d(samp,r=c(4-i,i),psi0=psi0)}
    Psi<-c(psi4s[1],2*psi4s[2],psi4s[3],4*psi4s[3],2*psi4s[4],psi4s[5])
    Psi<-invvech(Psi)
    return(Psi)
    }
    
###### Hpi selector (FULL)
    
pi.F<-function(samp,H){  ### The function PI(H)
    n<-nrow(samp)
    d<-ncol(samp)
    RK<-(4*pi)^(-d/2)
    psi.mat<-Psi4(samp)
    pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
    return(drop(pi.temp))
    }
    
Hpi.F<-function(samp,Hstart){ ### This gives the selector Hpi WITHOUT pre-sphering and PI(Hpi)
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    RK<-(4*pi)^(-d/2)
    if (missing(Hstart)) 
        Hstart <- matrix.sqrt((4/(n * (d + 2)))^(2/(d + 4)) * S)
    psi.mat<-Psi4(samp)
    pi.temp <- function(vechH) {
        H <- invvech(vechH) %*% invvech(vechH)
        pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% 
            psi.mat %*% vech(H)
        return(drop(pi.temp))
    }
    if (prod(eigen(psi.mat)$val > 0) == 1) {
        result <- optim(vech(Hstart), pi.temp, method = "BFGS")
        H <- invvech(result$par) %*% invvech(result$par)
        piH <- result$value
    }
    else {
        cat("Psi matrix not positive definite\n")
        H <- matrix(NA, nc = d, nr = d)
        piH<- NA
    }
    res<-list(H,piH)
    names(res)<-c("H","piH")
    return(res)
    }
    
Hpi.F2<-function(samp,Hstart){ ### This gives the selector Hpi WITH pre-sphering
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    S12<-matrix.sqrt(S)
    samp<-pre.sphere(samp)
    RK<-(4*pi)^(-d/2)
    if (missing(Hstart)) 
        Hstart <- sqrt((4/(n * (d + 2)))^(2/(d + 4))) * diag(d)
    psi.mat<-Psi4(samp)
    pi.temp <- function(vechH) {
        H <- invvech(vechH) %*% invvech(vechH)
        pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% 
            psi.mat %*% vech(H)
        return(drop(pi.temp))
    }
    if (prod(eigen(psi.mat)$val > 0) == 1) {
        result <- optim(vech(Hstart), pi.temp, method = "BFGS")
        H <- invvech(result$par) %*% invvech(result$par)
    }
    else {
        cat("Psi matrix not positive definite\n")
        H <- matrix(NA, nc = d, nr = d)
    }
    res<-S12 %*% H %*% S12
    return(res)
    }

###### Hpi selector (DIAGONAL)

Psi.D<-function(samp){
    Psi4(samp)
    }
    
pi.D<-function(samp,H){
    n<-nrow(samp)
    d<-ncol(samp)
    RK<-(4*pi)^(-d/2)
    psi.mat<-Psi.D(samp)
    H <- diag(H)
    pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
    return(drop(pi.temp))
    }
    
Hpi.D<-function(samp,Hstart){ ### Diagonal Hpi WITHOUT pre-scaling
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    RK<-(4*pi)^(-d/2)
    if (missing(Hstart)) 
        Hstart <- matrix.sqrt((4/(n * (d + 2)))^(2/(d + 4)) * S)
    psi.mat<-Psi.D(samp)
    pi.temp <- function(diagH) {
            H <- diag(diagH) %*% diag(diagH)
            pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
            return(drop(pi.temp))
        }
    result <- optim(diag(Hstart), pi.temp, method = "BFGS")
    H <- diag(result$par) %*% diag(result$par)
    piH <- result$value
    res<-list(H,piH)
    names(res)<-c("H","piH")
    return(res)
    }
        
Hpi.D2<-function(samp,Hstart){ ### Diagonal Hpi WITH pre-scaling
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    samp<-pre.scale(samp)
    RK<-(4*pi)^(-d/2)
    if (missing(Hstart)) 
        Hstart <- sqrt((4/(n * (d + 2)))^(2/(d + 4))) * diag(d)
    psi.mat<-Psi.D(samp)
    pi.temp <- function(diagH) {
            H <- diag(diagH) %*% diag(diagH)
            pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
            return(drop(pi.temp))
        }
    result <- optim(diag(Hstart), pi.temp, method = "BFGS")
    H <- diag(result$par) %*% diag(result$par)
    S12<-diag(sqrt(diag(S)))
    res<-S12 %*% H %*% S12
    return(res)
    }


###### Hpi selector (SINGLE)

Psi.S<-function(samp){
    Psi4(samp)
    }
    
pi.S<-function(samp,H){
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    RK<-(4*pi)^(-d/2)
    psi.mat<-Psi.D(samp)
    pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
    return(drop(pi.temp))
    }
    
Hpi.S<-function(samp,Hstart){ ### Single bandwith Hpi WITHOUT pre-sphering or pre-scaling
    n<-nrow(samp)
    d<-ncol(samp)
    S<-var(samp)
    RK<-(4*pi)^(-d/2)
    if (missing(Hstart)) 
        Hstart <- (4/(n * (d + 2)))^(2/(d + 4)) * S
    psi.mat<-Psi.D(samp)
    pi.temp <- function(h) {
            H <- h^2*diag(d)
            pi.temp <- 1/(det(H)^(1/2) * n) * RK + 1/4 * t(vech(H)) %*% psi.mat %*% vech(H)
            return(drop(pi.temp))
        }
    result <- optim(mean(diag(Hstart)), pi.temp, method = "BFGS")
    H <- (result$par)*diag(d)
    piH <- result$value
    res<-list(H,piH)
    names(res)<-c("H","piH")
    return(res)
    }  

#### These two functions give exactly the same selectors. Compare:
#### x<-rmvnorm(100,mean=rep(0,2),sigma=diag(2))
#### Hpi.F2(x) #My implementation
#### Hpi(x,pilot="a") #Tarn's implementation

#### Compare the mean implementation time along 20 samples

compare<-function(n){
    t1<-numeric(20)
    t2<-numeric(20)
    for(i in 1:20){
        x<-rmvnorm(n,mean=rep(0,2),sigma=diag(2))
        t1[i]<-system.time(Hpi.F2(x))[3]
        t2[i]<-system.time(Hpi(x,pilot="a"))[3]
        }
    return(c(mean(t1),mean(t2)))
    }
    
#### For n=100, Hpi(x,pilot="a") is more than 4 times slower than Hpi.F2(x)   
