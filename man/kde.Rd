\name{kde}
\alias{kda.kde}
\alias{kde}

\title{Kernel density estimate for multivariate data}
\description{
  Kernel density estimate for 1- to 6-dimensional data.
}
\usage{
kde(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
    binned=FALSE, bgridsize, positive=FALSE, adj.positive, w,
    compute.cont=FALSE, approx.cont=TRUE)

kda.kde(x, x.group, Hs, hs, prior.prob=NULL, gridsize, xmin, xmax,
    supp=3.7, eval.points=NULL, binned=FALSE, bgridsize, w,
    compute.cont=FALSE, approx.cont=TRUE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{x.group}{vector of group labels}
  \item{H,Hs}{bandwidth matrix(ces)}
  \item{h,hs}{scalar bandwidth(s)}
  \item{prior.prob}{vector of prior probabilities}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin}{vector of minimum values for grid}
  \item{xmax}{vector of maximum values for grid}
  \item{supp}{effective support for standard normal.}
  \item{eval.points}{points at which density estimate is evaluated}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
  \item{adj.positive}{adjustment applied to data \code{x <- log(x +
	adj.positive)} when \code{positive=TRUE}. Default is the minimum of \code{x}.}
  \item{w}{vector of weights (non-negative and sum is equal to sample size). Default is a vector of all ones.}
  \item{compute.cont}{flag for computing probability contour levels from 1\% to 99\%. Default is FALSE.}
  \item{approx.cont}{flag for computing approximate probability contour levels. Default is TRUE.}
}
  

\value{
  --The result from \code{kde} is a kernel density estimate which is an object of class \code{kde}: 
  \item{x}{data points - same as input}
  \item{eval.points}{points at which the density estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth (1-d only)}
  \item{w}{weights}
  \item{cont}{probability contour levels}

  --The result from \code{kda.kde} is a density estimate
  for discriminant analysis which is an object of class \code{kda.kde}:
  \item{x}{data points - same as input}
  \item{x.group}{group labels - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}  
  \item{prior.prob}{prior probabilities}
  \item{H}{bandwidth matrices (>1-d only) or} 
  \item{h}{bandwidths (1-d only)}
  \item{cont}{probability contour levels}
  \item{w}{weights}
}

\details{
  For d = 1, 2, 3, 4, and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid 
  defined by \code{gridsize} (if \code{binned=FALSE}) or
  by \code{bgridsize} (if \code{binned=TRUE}).
  If \code{eval.points} is specified, then the
  density estimate is computed exactly at \code{eval.points}.
  
  For d > 4, the kernel density estimate is computed exactly
  and \code{eval.points} must be specified.

  \code{supp} is the effective support for a normal kernel, i.e.  
  all values outside \code{[-supp,supp]^d} are set to zero. 
  The default \code{xmin} is \code{min(x)-Hmax*supp} and \code{xmax}
  is \code{max(x)+Hmax*supp}  where \code{Hmax} is the maximum of the
  diagonal elements of \code{H}.The default weights \code{w} is a vector of all ones.

  For \code{kda.kde}, if you have prior probabilities then set \code{prior.prob} to these.
  Otherwise the default is \code{prior.prob=NULL} i.e. use the sample
  proportions as estimates of the prior probabilities.
}

\seealso{\code{\link{plot.kde}}, \code{\link{plot.kda.kde}}}

\examples{
### See examples in ? plot.kde, ? plot.kda.kde  
}
\keyword{smooth}
