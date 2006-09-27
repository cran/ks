\name{kde}
\alias{kde}

\title{Kernel density estimate for multivariate data}
\description{
  Kernel density estimate for 1- to 6-dimensional data.
}
\usage{
 kde(x, H, h, gridsize, supp=3.7, eval.points)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth}
  \item{gridsize}{vector of number of grid points}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points that density estimate is evaluated at} 
}
  

\value{
  Kernel density estimate is an object of class \code{kde} which is a
  list with 4 fields
  \item{x}{data points - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrix (>1-d only) or}
  \item{h}{scalar bandwidth (1-d only)}
}

\details{
  For d > 1, the kernel density estimate is computed exactly i.e. binning
  is not used. For d = 1, the binned estimator from the
  \code{KernSmooth} library is used.

  For d = 1, 2, 3, if \code{eval.points} is not specified, then the
  density estimate is automatically computed over a grid whose
  resolution is controlled by \code{gridsize} (default is 101, 51 x 51
  and 51 x 51 x 51 respectively).

  For d > 3, \code{eval.points} must be specified.
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\seealso{\code{\link{plot.kde}}}

\examples{
### univariate example
data(unicef)
fhat <- kde(unicef[,1], h=sqrt(944))
fhat <- kde(unicef[,1], H=944) ## same as above

### bivariate example
data(unicef)
H.pi <- Hpi(unicef, nstage=1)
fhat <- kde(unicef, H=H.pi)

### 4-variate example
library(MASS)
data(iris)
ir <- iris[,1:4][iris[,5]=="setosa",]
H.scv <- Hscv(ir)
fhat <- kde(ir, H=H.scv, eval.points=ir)  
}
\keyword{smooth}
