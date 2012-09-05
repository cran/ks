\name{kda.kde}
\alias{kda.kde}

\title{Kernel density estimate for discriminant analysis for multivariate data}
\description{
  Kernel density estimate for discriminant analysis for 1- to 6-dimensional data.
}
\usage{
kda.kde(x, x.group, Hs, hs, prior.prob=NULL, gridsize, xmin, xmax,
    supp=3.7, eval.points=NULL, binned=FALSE, bgridsize, w,
    compute.cont=FALSE, approx.cont=TRUE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{x.group}{vector of group labels}
  \item{Hs}{bandwidth matrix(ces)}
  \item{hs}{scalar bandwidth(s)}
  \item{prior.prob}{vector of prior probabilities}
  \item{gridsize}{vector of number of grid points}
  \item{xmin}{vector of minimum values for grid}
  \item{xmax}{vector of maximum values for grid}
  \item{supp}{effective support for standard normal.}
  \item{eval.points}{points at which density estimate is evaluated}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights (non-negative and sum is equal to sample size). Default is a vector of all ones.}
  \item{compute.cont}{flag for computing probability contour levels from 1\% to 99\%. Default is FALSE.}
  \item{approx.cont}{flag for computing approximate probability contour levels. Default is TRUE.}
}
  

\value{
  A density estimate for discriminant analysis is an object of class \code{kda.kde}:
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
  If you have prior probabilities then set \code{prior.prob} to these.
  Otherwise the default is \code{prior.prob=NULL} i.e. use the sample
  proportions as estimates of the prior probabilities.
  
  For details of the other parameters, see \code{\link{kde}}.
}

\seealso{\code{\link{plot.kda.kde}}}

\examples{
## trivariate example
## colour indicates species, transparency indicates density heights
library(MASS)
data(iris)
ir <- iris[,1:3]
ir.gr <- iris[,5] 
H <- Hkda(ir, ir.gr, bw="plugin", pilot="dscalar")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H)
plot(kda.fhat)
}
\keyword{smooth}
