

.First.lib <- function(lib=NULL, pkg=ks)
{
  library.dynam("ks", pkg, lib)
  x <- installed.packages()
  cat(pkg, x[x[,1]==pkg,3], "(2006)\n")
  rm(x)
  require(mvtnorm)
  require(misc3d)
  require(rgl)
  #require(feature)
}  

