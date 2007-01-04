\name{kde}
\alias{kde}

\title{Kernel density estimate for multivariate data}
\description{
  Kernel density estimate for 1- to 6-dimensional data.
}
\usage{
 kde(x, H, h, gridsize, binned=FALSE, bgridsize, supp=3.7, eval.points)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth}
  \item{gridsize}{vector of number of grid points - used for display and
    binning}
  \item{binned}{flag for binned estimation (default is FALSE)}
  \item{bgridsize}{vector of binning grid sizes - only required if \code{binned=TRUE}}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points at which density estimate is evaluated} 
}
  

\value{
  Kernel density estimate is an object of class \code{kde} which is a
  list with 4 fields
  \item{x}{data points - same as input}
  \item{eval.points}{points at which the density estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrix (>1-d only) or}
  \item{h}{scalar bandwidth (1-d only)}
}

\details{
  For d = 1, 2, 3, 4, 
  and if \code{eval.points} is not specified, then the
  density estimate is computed over a grid 
  defined by \code{gridsize} (if \code{binned=FALSE}) or
  by \code{bgridsize} (if \code{binned=TRUE}).

  For d = 1, 2, 3, 4, 
  and if \code{eval.points} is specified, then the
  density estimate is computed exactly at \code{eval.points}.
  
  For d > 4, the kernel density estimate is computed exactly
  and \code{eval.points} must be specified.
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\seealso{\code{\link{plot.kde}}}

\examples{
### See examples in ? plot.kde  
}
\keyword{smooth}
