\name{kda.kde}
\alias{kda.kde}

\title{Kernel density estimate for kernel discriminant analysis for multivariate data}
\description{
  Kernel density estimate for kernel discriminant analysis for 1- to 6-dimensional data
}
\usage{
kda.kde(x, x.group, Hs, hs, prior.prob=NULL, gridsize, supp=3.7,
        eval.points=NULL)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{Hs}{(stacked) matrix of bandwidth matrices}
  \item{hs}{vector of scalar bandwidths}
  \item{prior.prob}{vector of prior probabilities}
  \item{gridsize}{vector of number of grid points}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points that density estimate is evaluated at}
}
  

\value{
  The kernel density estimate for kernel discriminant analysis is
  based on \code{\link{kde}}, one density estimate for each group. 

  The result from \code{kda.kde} is a density estimate
  for discriminant analysis is an object of class \code{kda.kde} which is a
  list with 6 fields
  \item{x}{data points - same as input}
  \item{x.group}{group labels - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}  
  \item{prior.prob}{sample proportions of each group}
  \item{H}{bandwidth matrices (>1-d only)  or } 
  \item{h}{bandwidths (1-d only)}
}

\details{
  If you have prior probabilities then set \code{prior.prob} to these.
  Otherwise \code{prior.prob=NULL} is the default i.e. use the sample
  proportions as estimates of the prior probabilities.

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

\seealso{\code{\link{plot.kda.kde}}}

\examples{
### bivariate example - restricted iris dataset  
library(MASS)
data(iris)
ir <- iris[,1:2]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.fhat <- kda.kde(ir, ir.gr, H=H)
}
\keyword{smooth}
