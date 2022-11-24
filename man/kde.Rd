\name{kde}
\alias{kde}
\alias{predict.kde}
\title{Kernel density estimate}

\description{
  Kernel density estimate for 1- to 6-dimensional data.
}

\usage{
kde(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, 
    bgridsize, positive=FALSE, adj.positive, w, compute.cont=TRUE, 
    approx.cont=TRUE, unit.interval=FALSE, density=FALSE, verbose=FALSE)

\method{predict}{kde}(object, ..., x, zero.flag=TRUE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{vector or matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if data are positive (1-d, 2-d). Default is FALSE.}
  \item{adj.positive}{adjustment applied to positive 1-d data}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{compute.cont}{flag for computing 1\% to 99\% probability contour levels. Default is TRUE.}
  \item{approx.cont}{flag for computing approximate probability contour
    levels. Default is TRUE.}
  \item{unit.interval}{flag if 1-d data are bounded on unit interval
  [0,1]. Default is FALSE.} 
  \item{density}{flag if density estimate values are forced to be non-negative function. Default is FALSE.} 
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{object}{object of class \code{kde}}
  \item{zero.flag}{deprecated (retained for backwards compatibilty)}
  \item{...}{other parameters}
}

\value{
  A kernel density estimate is an object of class \code{kde} which is a
  list with fields: 
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{w}{vector of weights}
  \item{cont}{vector of probability contour levels}    
}

\details{
  For d=1, if \code{h} is missing, the default bandwidth is \code{hpi}.
  For d>1, if \code{H} is missing, the default is \code{Hpi}.

  For d=1, if \code{positive=TRUE} then \code{x} is transformed to
  \code{log(x+adj.positive)} where the default \code{adj.positive} is
  the minimum of \code{x}. This is known as a log transformation density
  estimate. 
  
  For d=1, 2, 3, and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid 
  defined by \code{gridsize} (if \code{binned=FALSE}) or
  by \code{bgridsize} (if \code{binned=TRUE}). This form is suitable for
  visualisation in conjunction with the \code{plot} method.

  For d=4, 5, 6,  and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid defined by \code{gridsize}. 
  
  If \code{eval.points} is specified, as a vector (d=1) or 
  as a matrix (d=2, 3, 4), then the density estimate is computed at
  \code{eval.points}. This form is suitable for numerical summaries
  (e.g. maximum likelihood), and is not compatible with the \code{plot}
  method. Despite that the density estimate is returned only at
  \code{eval.points}, by default, a binned gridded estimate is
  calculated first and then the density estimate at \code{eval.points}
  is computed using the \code{predict} method. If this default intermediate
  binned grid estimate is not required, then set \code{binned=FALSE} to
  compute directly the exact density estimate at \code{eval.points}.   

  Binned kernel estimation is an approximation to the exact kernel
  estimation and is available for d=1, 2, 3, 4. This makes
  kernel estimators feasible for large samples. The default value of the
  binning flag \code{binned} is n>1 (d=1), n>500 (d=2), n>1000 (d>=3). 
  Some times binned estimation leads to negative density values: if non-negative
  values are required, then set \code{density=TRUE}. 
  
  The default \code{bgridsize,gridsize} are d=1: 401; d=2: rep(151, 2); d=3: 
  rep(51, 3); d=4: rep(21, 4). 
  
  The effective support for a normal kernel is where  
  all values outside \code{[-supp,supp]^d} are zero. 

  The default \code{xmin} is \code{min(x)-Hmax*supp} and \code{xmax}
  is \code{max(x)+Hmax*supp}  where \code{Hmax} is the maximum of the
  diagonal elements of \code{H}. The grid produced is the outer
  product of \code{c(xmin[1], xmax[1])}, ..., \code{c(xmin[d], xmax[d])}. 
  For \pkg{ks} \eqn{\geq}{>=} 1.14.0, when \code{binned=TRUE} and \code{xmin,xmax} 
  are not missing, the data values \code{x} are clipped to the estimation grid 
  delimited by \code{xmin,xmax} to prevent potential memory leaks.
}

\seealso{\code{\link{plot.kde}}, \code{\link{kde.boundary}}}

\examples{
## unit interval data 
set.seed(8192)             
fhat <- kde(runif(10000,0,1), unit.interval=TRUE)
plot(fhat, ylim=c(0,1.2))

## positive data 
data(worldbank)
wb <- as.matrix(na.omit(worldbank[,2:3]))
wb[,2] <- wb[,2]/1000
fhat <- kde(x=wb)
fhat.trans <- kde(x=wb, adj.positive=c(0,0), positive=TRUE)
plot(fhat, col=1, xlim=c(0,20), ylim=c(0,80))
plot(fhat.trans, add=TRUE, col=2)
rect(0,0,100,100, lty=2)

## large data on non-default grid
## 151 x 151 grid = [-5,-4.933,..,5] x [-5,-4.933,..,5]
set.seed(8192)
x <- rmvnorm.mixt(10000, mus=c(0,0), Sigmas=invvech(c(1,0.8,1)))
fhat <- kde(x=x, compute.cont=TRUE, xmin=c(-5,-5), xmax=c(5,5), bgridsize=c(151,151))
plot(fhat)

## See other examples in ? plot.kde
}

\keyword{smooth}
