\name{kfs}
\alias{kfs}
\title{Kernel feature significance }

\description{
  Kernel feature significance for 1- to 6-dimensional data.
}

\usage{
kfs(x, H, h, deriv.order=2, gridsize, gridtype, xmin, xmax, supp=3.7,
    eval.points, binned, bgridsize, positive=FALSE, adj.positive, w, 
    verbose=FALSE, signif.level=0.05)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{deriv.order}{derivative order (scalar)}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{vector or matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
  \item{adj.positive}{adjustment applied to positive 1-d data}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{signif.level}{overall level of significance for hypothesis
  tests. Default is 0.05.} 
}
  
\value{
  A kernel feature significance estimate is an object of class
  \code{kfs} which is a list with fields 
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{binary matrix for significant feature at
  \code{eval.points}: 0 = not signif.,  1 = signif.}
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names}
  \item{w}{vector of weights}
  \item{deriv.order}{derivative order (scalar)}
  \item{deriv.ind}{martix where each row is a vector of partial derivative indices.}

  This is the same structure as a \code{kdde} object, except that
  \code{estimate} is a binary matrix rather than real-valued. 
}

\details{
  Feature significance is based on significance testing of the gradient
  (first derivative) and curvature (second derivative) of a kernel
  density estimate. Only the latter is currently implemented, and is
  also known as significant modal regions.

  The hypothesis test at a grid point \eqn{\bold{x}}{x} is
  \eqn{H_0(\bold{x}): \mathsf{H} f(\bold{x}) < 0}{H0(x): H f(x) < 0},
  i.e.  the density Hessian matrix \eqn{\mathsf{H} f(\bold{x})}{H f(x)} is negative definite.
  The \eqn{p}{p}-values are computed for each \eqn{\bold{x}}{x} using that
  the test statistic is
  approximately chi-squared distributed with \eqn{d(d+1)/2}{d(d+1)/2} d.f.
  We then use a Hochberg-type simultaneous testing procedure, based on the
  ordered \eqn{p}{p}-values, to control the
  overall level of significance to be \code{signif.level}. If
  \eqn{H_0(\bold{x})}{H0(x)} is rejected then \eqn{\bold{x}}{x}
  belongs to a significant modal region. 

  The computations are based on \code{kdde(x, deriv.order=2)} so
  \code{kfs} inherits its behaviour from \code{\link{kdde}}.
  If the bandwidth \code{H} is missing, then
  the default bandwidth is the plug-in selector
  \code{Hpi(,deriv.order=2)}. Likewise for missing \code{h}.
  The effective support, binning, grid size, grid range, positive
  parameters are the same as \code{\link{kde}}.

  This function is similar to the \code{featureSignif} function in the
  \pkg{feature} package, except that it accepts unconstrained bandwidth
  matrices. 
}

\references{
  Chaudhuri, P. & Marron, J.S. (1999) 
  SiZer for exploration of structures in curves.
  \emph{Journal of the American Statistical Association},
  \bold{94},  807-823.

  Duong, T., Cowling, A., Koch, I. & Wand, M.P. (2008)
  Feature significance for multivariate kernel density estimation.
  \emph{Computational Statistics and Data Analysis}, \bold{52},
  4225-4242. 
  
  Godtliebsen, F., Marron, J.S. & Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.
}

\seealso{\code{\link{kdde}}, \code{\link{plot.kfs}}}

\examples{
library(MASS)
data(geyser)
geyser.fs <- kfs(geyser$duration, binned=TRUE)
plot(geyser.fs, xlab="duration")

## see example in ? plot.kfs
}

\keyword{smooth}
