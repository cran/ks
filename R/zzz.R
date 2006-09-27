

##.First.lib <- function(lib=NULL, pkg=ks)
#{
#  library.dynam("ks", pkg, lib)
  ##x <- installed.packages()
  ##cat(pkg, x[x[,1]==pkg,3], "(2006)\n")
  ##rm(x)
#  require(mvtnorm)
#  require(misc3d)
#  require(rgl)
#}  

.onLoad <- function(libname=NULL, pkgname=ks)
{
  #x <- installed.packages()
  #cat(pkgname, x[x[,1]==pkgname,3], "(2006)\n")
  #rm(x)
  cat("ks 1.4.3 (2006)\n")
}


