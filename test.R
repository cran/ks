


for (ff in system("ls R/*.R", intern=TRUE)) source(ff)

mu1 <- rep(0,5)
mu2 <- rep(1,5,5)
Sigma1 <- 0.16*diag(c(1,2, 1,2,1))
Sigma2 <- 0.16*diag(5)
props <- c(1,1)/2

set.seed(8192)
y <- rmvnorm.mixt(n=1000, mus=rbind(mu1, mu2), Sigmas=rbind(Sigma1, Sigma2), props=props)

dmvnorm.mixt(x=, mus=rbind(mu1, mu2), Sigmas=rbind(Sigma1, Sigma2), props=props)

#unix.time(H <- Hpi(y[,1:2], nstage=1, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:2], nstage=1, binned=TRUE)); print(H)
#unix.time(H <- Hpi(y[,1:3], nstage=1, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:3], nstage=1, binned=TRUE)); print(H)
#unix.time(H <- Hpi(y[,1:4], nstage=1, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:4], nstage=1, binned=TRUE)); print(H)

#unix.time(H <- Hpi(y[,1:2], nstage=2, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:2], nstage=2, binned=TRUE)); print(H)
#unix.time(H <- Hpi(y[,1:3], nstage=2, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:3], nstage=2, binned=TRUE)); print(H)
#unix.time(H <- Hpi(y[,1:4], nstage=2, binned=FALSE)); print(H)
#unix.time(H <- Hpi(y[,1:4], nstage=2, binned=TRUE)); print(H)

#unix.time(H <- Hpi.diag(y[,1:2], nstage=1, binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:2], nstage=1, binned=TRUE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:3], nstage=1, pilot="samse", binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:3], nstage=1, pilot="samse", binned=TRUE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:4], nstage=1, pilot="samse", binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:4], nstage=1, pilot="samse", binned=TRUE)); print(H)

#unix.time(H <- Hpi.diag(y[,1:2], nstage=2, binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:2], nstage=2, binned=TRUE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:3], nstage=2, pilot="samse", binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:3], nstage=2, pilot="samse", binned=TRUE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:4], nstage=2, pilot="samse", binned=FALSE)); print(H)
#unix.time(H <- Hpi.diag(y[,1:4], nstage=2, pilot="samse", binned=TRUE)); print(H)


#unix.time(H <- Hscv(y[,1:2], binned=FALSE)); print(H)
#unix.time(H <- Hscv(y[,1:2], binned=TRUE)); print(H)
#unix.time(H <- Hscv(y[,1:3], binned=FALSE)); print(H)
#unix.time(H <- Hscv(y[,1:3], binned=TRUE)); print(H)
#unix.time(H <- Hscv(y[,1:4], binned=FALSE)); print(H)
#unix.time(H <- Hscv(y[,1:4], binned=TRUE)); print(H)


unix.time(H <- Hlscv(y[,1:2])); print(H)
unix.time(H <- Hlscv(y[,1:3])); print(H)
unix.time(H <- Hlscv(y[,1:4])); print(H)


unix.time(H <- Hlscv.diag(y[,1:2], binned=FALSE)); print(H)
unix.time(H <- Hlscv.diag(y[,1:2], binned=TRUE)); print(H)
unix.time(H <- Hlscv.diag(y[,1:3], binned=FALSE)); print(H)
unix.time(H <- Hlscv.diag(y[,1:3], binned=TRUE)); print(H)
unix.time(H <- Hlscv.diag(y[,1:4], binned=FALSE)); print(H)
unix.time(H <- Hlscv.diag(y[,1:4], binned=TRUE)); print(H)
