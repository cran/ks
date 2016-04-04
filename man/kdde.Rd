\name{kdde}
\alias{kdde}
\alias{predict.kdde}

\title{Kernel density derivative estimate}
\description{
  Kernel density derivative estimate for 1- to 6-dimensional data.
}
\usage{
kdde(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, 
    eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, w,
    deriv.vec=TRUE, verbose=FALSE)

\method{predict}{kdde}(object, ..., x)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{deriv.order}{derivative order (scalar)}
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
  \item{deriv.vec}{flag to compute all derivatives in vectorised
  derivative. Default is TRUE. If FALSE then only the unique derivatives are computed.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.} 
  \item{object}{object of class \code{kdde}}
  \item{...}{other parameters} 
}
  

\value{
  A kernel density derivative estimate is an object of class
  \code{kdde} which is a list with fields:
  \item{x}{data points - same as input}
  \item{eval.points}{points at which the estimate is evaluated}
  \item{estimate}{density derivative estimate at \code{eval.points}}
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names}
  \item{w}{weights}
  \item{deriv.order}{derivative order (scalar)}
  \item{deriv.ind}{each row is a vector of partial derivative indices}
}

\details{ 
  For each partial derivative, for grid estimation, the estimate is a
  list whose elements  
  correspond to the partial derivative indices in the rows of \code{deriv.ind}. 
  For points estimation, the estimate is a matrix whose columns correspond to 
  rows of \code{deriv.ind}.

  If the bandwidth \code{H} is missing from \code{kdde}, then
  the default bandwidth is the  plug-in selector
  \code{Hpi}. Likewise for missing \code{h}.
  
  The effective support, binning, grid size, grid range, positive
  parameters are the same as \code{\link{kde}}. 
}

\seealso{\code{\link{kde}}}

\examples{
set.seed(8192)
x <- rmvnorm.mixt(1000, mus=c(0,0), Sigmas=invvech(c(1,0.8,1)))
fhat <- kdde(x=x, binned=TRUE, deriv.order=1) ## gradient [df/dx, df/dy]
predict(fhat, x=x[1:5,])

## See other examples in ? plot.kdde
}

\keyword{smooth}
