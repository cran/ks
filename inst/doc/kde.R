### R code from vignette source 'kde.Rnw'

###################################################
### code chunk number 1: kde.Rnw:99-107
###################################################
library(ks)
set.seed(8192)
samp <- 200
mus <- rbind(c(-2,2), c(0,0), c(2,-2))
Sigmas <- rbind(diag(2), matrix(c(0.8, -0.72, -0.72, 0.8), nrow=2), diag(2))
cwt <- 3/11
props <- c((1-cwt)/2, cwt, (1-cwt)/2)
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)


###################################################
### code chunk number 2: kde.Rnw:113-114
###################################################
plotmixt(mus=mus, Sigmas=Sigmas, props=props, xlim=c(-4,4), ylim=c(-4,4))


###################################################
### code chunk number 3: kde.Rnw:116-117
###################################################
plot(x, xlim=c(-4,4), ylim=c(-4,4), xlab="x", ylab="y")


###################################################
### code chunk number 4: kde.Rnw:126-128
###################################################
Hpi1 <- Hpi(x=x)
Hpi2 <- Hpi.diag(x=x)


###################################################
### code chunk number 5: kde.Rnw:132-134
###################################################
fhat.pi1 <- kde(x=x, H=Hpi1)
fhat.pi2 <- kde(x=x, H=Hpi2)


###################################################
### code chunk number 6: kde.Rnw:142-144 (eval = FALSE)
###################################################
## plot(fhat.pi1)
## plot(fhat.pi2)


###################################################
### code chunk number 7: kde.Rnw:156-157
###################################################
plot(fhat.pi1, main="Plug-in", cex.main=1.4, xlim=c(-4,4), ylim=c(-4,4))


###################################################
### code chunk number 8: kde.Rnw:159-160
###################################################
plot(fhat.pi2, main="Plug-in diagonal", cex.main=1.4, xlim=c(-4,4), ylim=c(-4,4)) 


###################################################
### code chunk number 9: kde.Rnw:171-173
###################################################
Hscv1 <- Hscv(x=x)
Hscv2 <- Hscv.diag(x=x)


###################################################
### code chunk number 10: kde.Rnw:177-179
###################################################
fhat.cv1 <- kde(x=x, H=Hscv1)
fhat.cv2 <- kde(x=x, H=Hscv2)


###################################################
### code chunk number 11: kde.Rnw:181-182
###################################################
plot(fhat.cv1, main="SCV", cex.main=1.4, xlim=c(-4,4), ylim=c(-4,4))


###################################################
### code chunk number 12: kde.Rnw:184-185
###################################################
plot(fhat.cv2, main="SCV diagonal", cex.main=1.4, xlim=c(-4,4), ylim=c(-4,4))


