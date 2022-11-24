\name{kcopula}
\alias{kcopula}
\alias{kcopula.de}
\title{Kernel copula (density) estimate}

\description{
  Kernel copula and copula density estimator for 2-dimensional data.
}

\usage{
kcopula(x, H, hs, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points,
  binned, bgridsize, w, marginal="kernel", verbose=FALSE)
kcopula.de(x, H, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
  binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE,
  marginal="kernel", boundary.supp, boundary.kernel="beta", verbose=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,hs}{bandwidth matrix. If these are missing, \code{Hpi.kcde}/\code{Hpi} or \code{hpi.kcde}/\code{hpi} is called by default.}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{marginal}{"kernel" = kernel cdf or "empirical" = empirical cdf
    to calculate pseudo-uniform values. Default is "kernel".}
  \item{compute.cont}{flag for computing 1\% to 99\% probability contour levels. Default is TRUE.}
  \item{approx.cont}{flag for computing approximate probability contour
    levels. Default is TRUE.}
  \item{boundary.supp}{effective support for boundary region}
  \item{boundary.kernel}{"beta" = beta boundary kernel, "linear" = linear boundary kernel}
 \item{verbose}{flag to print out progress information. Default is FALSE.}  
}
  
\value{
  A kernel copula estimate, output from \code{kcopula}, is an object of
  class \code{kcopula}. A kernel copula density estimate, output from
  \code{kcopula.de}, is an object of class \code{kde}. These two classes
  of objects have the same fields as \code{kcde} and \code{kde} objects
  respectively, except for
  \item{x}{pseudo-uniform data points}
  \item{x.orig}{data points - same as input}
  \item{marginal}{marginal function used to compute pseudo-uniform data}
  \item{boundary}{flag for data points in the boundary region
  (\code{kcopula.de} only)}
}

\details{
  For kernel copula estimates, a transformation approach is used to
  account for the boundary effects. If \code{H} is missing, the default
  is \code{Hpi.kcde}; if \code{hs} are missing, the default is
  \code{hpi.kcde}.  
  
  For kernel copula density estimates, for those points which are in
  the interior region, the usual kernel density estimator
  (\code{\link{kde}}) is used. For those points in the boundary region,
  a product beta kernel based on the boundary corrected univariate beta
  kernel of Chen (1999) is used (\code{\link{kde.boundary}}). If \code{H}
  is missing, the default is \code{Hpi.kcde}; if \code{hs} are missing,
  the default is \code{hpi}.
  
  The effective support, binning, grid size, grid range parameters are
  the same  as for \code{\link{kde}}.  
  
}
 
\references{
  Duong, T. (2014) Optimal data-based smoothing for non-parametric
  estimation of copula functions and their densities. Submitted.

  Chen, S.X. (1999). Beta kernel estimator for density
  functions. \emph{Computational Statistics & Data Analysis},
  \bold{31}, 131--145.
}

\seealso{\code{\link{kcde}}, \code{\link{kde}}}

\examples{
library(MASS)
data(fgl)
x <- fgl[,c("RI", "Na")]
Chat <- kcopula(x=x)
plot(Chat, display="persp", border=1)
plot(Chat, display="filled.contour", lwd=1)
}

\keyword{smooth}
