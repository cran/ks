\name{kdde}
\alias{kdde}

\title{Kernel density derivative estimate for multivariate data}
\description{
  Kernel density derivative estimate for 1- to 6-dimensional data.
}
\usage{
kdde(x, H, h, deriv.order=0, gridsize, gridtype, xmin, xmax, supp=3.7, 
    eval.points, binned=FALSE, bgridsize, positive=FALSE, adj.positive, 
    w, deriv.vec=TRUE, verbose=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth}
  \item{deriv.order}{derivative order (scalar)}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin}{vector of minimum values for grid}
  \item{xmax}{vector of maximum values for grid}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points at which density estimate is evaluated}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
  \item{adj.positive}{adjustment applied to data \code{x <- log(x +
	adj.positive)} when \code{positive=TRUE}. Default is the minimum of \code{x}.}
  \item{w}{vector of weights (non-negative and sum is equal to sample size). Default is a vector of all ones.}
  \item{deriv.vec}{flag to compute all derivatives in vectorised derivative. Default is TRUE. 
       If FALSE then only the unique derivatives are computed.}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
}
  

\value{
  The result is kernel density derivative estimate is an object of class \code{kdde}:
  \item{x}{data points - same as input}
  \item{eval.points}{points at which the density estimate is evaluated}
  \item{estimate}{density derivative estimate at \code{eval.points}}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth (1-d only)}
  \item{w}{weights}
  \item{deriv.order}{derivative order (scalar)}
  \item{deriv.ind}{each row is a vector of partial derivative indices}
}

\details{ 
  Binning/exact estimation behaviour is the same as for \code{\link{kde}}.

  For each partial derivative, for grid estimation, the estimate is a list whose elements 
  correspond to the partial derivative indices in the rows of \code{deriv.ind}. 
  For points estimation, the estimate is a matrix whose columns correspond to 
  rows of \code{deriv.ind}. 
}

\seealso{\code{\link{kde}}}

\examples{
## univariate example
x <- rnorm.mixt(n=100, mus=1, sigmas=1, props=1)
fhat2 <- kdde(x=x, h=hpi(x), deriv.order=2)      ## d^2 f/dx^2 

## bivariate example
data(unicef)
H.scv <- Hscv(x=unicef)
fhat1 <- kdde(x=unicef, H=H.scv, deriv.order=1)   ## gradient vector
}
\keyword{smooth}
