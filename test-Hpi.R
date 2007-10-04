
for (ff in system("ls R/*.R", intern=TRUE)) source(ff);

source("HPI.R")
library(ks)
d <- 3
x <- rmvnorm.mixt(1000, mus=rep(0,d), diag(d))

H <- Hpi(x, binned=FALSE)