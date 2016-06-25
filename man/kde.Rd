\name{kde}
\alias{kde}
\alias{predict.kde}

\title{Kernel density estimate}

\description{
  Kernel density estimate for 1- to 6-dimensional data.
}
\usage{
kde(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
    binned=FALSE, bgridsize, positive=FALSE, adj.positive, w,
    compute.cont=FALSE, approx.cont=TRUE, unit.interval=FALSE,
    verbose=FALSE)

\method{predict}{kde}(object, ..., x, zero.flag=TRUE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{points at which estimate is evaluated}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
   \item{adj.positive}{adjustment applied to positive 1-d data}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{compute.cont}{flag for computing 1\% to 99\% probability contour levels. Default is FALSE.}
  \item{approx.cont}{flag for computing approximate probability contour
    levels. Default is TRUE.}
  \item{unit.interval}{flag if 1-d data are bounded on unit interval
  [0,1]. Default is FALSE.} 
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{object}{object of class \code{kde}}
  \item{zero.flag}{interploted values \code{object$estimate} of when \code{x}
  outside of \code{object$eval.points}=0 (if \code{TRUE}), = nearest
  \code{object$estimate} (if \code{FALSE})}
  \item{...}{other parameters}
}
  

\value{
  A kernel density estimate is an object of class \code{kde} which is a
  list with fields: 
  \item{x}{data points - same as input}
  \item{eval.points}{points at which the estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{w}{weights}
  \item{cont}{probability contour levels (if \code{compute.cont=TRUE})}    
}

\details{
  For d=1, if \code{h} is missing, the default bandwidth is \code{hpi}.
  For d>1, if \code{H} is missing, the default is \code{Hpi}.

  For d=1, if \code{positive=TRUE} then \code{x<-log(x+adj.positive)}
  where the default \code{adj.positive} is the minimum of \code{x}.
  
  For d=1, 2, 3, 4, and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid 
  defined by \code{gridsize} (if \code{binned=FALSE}) or
  by \code{bgridsize} (if \code{binned=TRUE}).
  If \code{eval.points} is specified, then the
  density estimate is computed exactly at \code{eval.points}.
  For d>4, the kernel density estimate is computed exactly
  and \code{eval.points} must be specified.

  The default \code{bgridsize,gridsize} are  d=1: 401; d=2: rep(151, 2);
  d=3: rep(51, 3); d=4: rep(21,4).
  
  The effective support for a normal kernel is where  
  all values outside \code{[-supp,supp]^d} are zero. 

  The default \code{xmin} is \code{min(x)-Hmax*supp} and \code{xmax}
  is \code{max(x)+Hmax*supp}  where \code{Hmax} is the maximum of the
  diagonal elements of \code{H}. The grid produced is the outer
  product of \code{c(xmin[1], xmax[1])}, ..., \code{c(xmin[d], xmax[d])}.
}

\seealso{\code{\link{plot.kde}}}

\examples{
## positive data example
set.seed(8192)
x <- 2^rnorm(100)
fhat <- kde(x=x, positive=TRUE)
plot(fhat, col=3)
points(c(0.5, 1), predict(fhat, x=c(0.5, 1)))

## large data example on non-default grid
## 151 x 151 grid = [-5,-4.933,..,5] x [-5,-4.933,..,5]
set.seed(8192)
x <- rmvnorm.mixt(10000, mus=c(0,0), Sigmas=invvech(c(1,0.8,1)))
fhat <- kde(x=x, binned=TRUE, compute.cont=TRUE, xmin=c(-5,-5), xmax=c(5,5), bgridsize=c(151,151))
plot(fhat)

## See other examples in ? plot.kde
}
\keyword{smooth}
