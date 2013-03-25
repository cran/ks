\name{kde}
\alias{kde}

\title{Kernel density estimate for multivariate data}
\description{
  Kernel density estimate for 1- to 6-dimensional data.
}
\usage{
kde(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
    binned=FALSE, bgridsize, positive=FALSE, adj.positive, w,
    compute.cont=FALSE, approx.cont=TRUE)
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
  \item{approx.cont}{flag for computing approximate probability contour levels. Default is TRUE.}
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
   For d=1, if \code{h} is missing, the default bandwidth is \code{hpi(x=x, nstage=2, binned=TRUE)},
  and for d>1, if \code{H} is missing, the default is \code{Hpi(x=x, nstage=2, binned=n>1000, pilot="samse")}.

  \code{supp} is the effective support for a normal kernel, i.e.  
  all values outside \code{[-supp,supp]^d} are set to zero. 
  The default \code{xmin} is \code{min(x)-Hmax*supp} and \code{xmax}
  is \code{max(x)+Hmax*supp}  where \code{Hmax} is the maximum of the
  diagonal elements of \code{H}.The default weights \code{w} is a vector
  of all ones.

  For d=1, 2, 3, 4, and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid 
  defined by \code{gridsize} (if \code{binned=FALSE}) or
  by \code{bgridsize} (if \code{binned=TRUE}).
  If \code{eval.points} is specified, then the
  density estimate is computed exactly at \code{eval.points}.
  For d>4, the kernel density estimate is computed exactly
  and \code{eval.points} must be specified.

  For d=1, if \code{positive=TRUE} then \code{x<-log(x+adj.positive)}
  where the default \code{adj.positive} is the minimum of \code{x}.
}

\seealso{\code{\link{plot.kde}}}

\examples{
## positive data example
x <- rnorm(100)
z <- 2^x
fhatz <- kde(x=z, positive=TRUE)
plot(fhatz)

## See other examples in ? plot.kde
}
\keyword{smooth}
