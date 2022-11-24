\name{kcde}
\alias{kcde}
\alias{Hpi.kcde}
\alias{Hpi.diag.kcde}
\alias{hpi.kcde}
\alias{predict.kcde}
\title{Kernel cumulative distribution/survival function estimate}

\description{
  Kernel cumulative distribution/survival function estimate for 1- to 3-dimensional data.
}

\usage{
kcde(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, binned, 
  bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE, 
  tail.flag="lower.tail")
Hpi.kcde(x, nstage=2, pilot, Hstart, binned, bgridsize, amise=FALSE, 
  verbose=FALSE, optim.fun="optim", pre=TRUE)
Hpi.diag.kcde(x, nstage=2, pilot, Hstart, binned, bgridsize, amise=FALSE,
  verbose=FALSE, optim.fun="optim", pre=TRUE)
hpi.kcde(x, nstage=2, binned, amise=FALSE)

\method{predict}{kcde}(object, ..., x)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, then
  \code{Hpi.kcde} or \code{hpi.kcde} is called by default.}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{vector or matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
  \item{adj.positive}{adjustment applied to positive 1-d data}
  \item{w}{not yet implemented}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{tail.flag}{"lower.tail" = cumulative distribution, "upper.tail" =
  survival function}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{"dscalar" = single pilot bandwidth (default for
  \code{Hpi.diag.kcde} \cr
    "dunconstr" = single unconstrained pilot bandwidth (default for
  \code{Hpi.kcde}} 
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{optim.fun}{optimiser function: one of \code{nlm} or
    \code{optim}}
  \item{pre}{flag for pre-scaling data. Default is TRUE.}
  \item{object}{object of class \code{kcde}}
  \item{...}{other parameters} 
}
  
\value{
  A kernel cumulative distribution estimate is an object of class
  \code{kcde} which is a list with fields:
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{cumulative distribution/survival function estimate at
  \code{eval.points}} 
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names}
  \item{w}{vector of weights}
  \item{tail}{"lower.tail"=cumulative distribution, "upper.tail"=survival function}
}

\details{
  If \code{tail.flag="lower.tail"} then the cumulative distribution
  function \eqn{\mathrm{Pr}(\bold{X}\leq\bold{x})}{Pr(X<=x)} is estimated, otherwise
  if \code{tail.flag="upper.tail"}, it is the survival function
  \eqn{\mathrm{Pr}(\bold{X}>\bold{x})}{P(X>x)}. For \eqn{d>1}{d>1},
  \eqn{\mathrm{Pr}(\bold{X}\leq\bold{x}) \neq 1 - \mathrm{Pr}(\bold{X}>\bold{x})}{Pr(X<=x) != 1-Pr(X>x)}.
  
  If the bandwidth \code{H} is missing in \code{kcde}, then
  the default bandwidth is the plug-in selector
  \code{Hpi.kcde}. Likewise for missing \code{h}.
  No pre-scaling/pre-sphering is used since the \code{Hpi.kcde} is not
  invariant to translation/dilation.

  The effective support, binning, grid size, grid range, positive, optimisation function
  parameters are the same as \code{\link{kde}}.
}

\references{
  Duong, T. (2016) Non-parametric smoothed estimation of multivariate
  cumulative distribution and survival functions, and receiver operating
  characteristic curves. \emph{Journal of the Korean Statistical
  Society}, \bold{45}, 33-50.
}

\seealso{\code{\link{kde}}, \code{\link{plot.kcde}}}

\examples{
data(iris, package="MASS")
Fhat <- kcde(iris[,1:2])  
predict(Fhat, x=as.matrix(iris[,1:2]))

## See other examples in ? plot.kcde
}

\keyword{smooth}
