\name{kroc}
\alias{kroc}
\alias{summary.kroc}

\title{Kernel receiver operating characteristic (ROC) curve}

\description{ Kernel receiver operating characteristic (ROC) curve for 1- to 3-dimensional data.}

\usage{
kroc(x1, x2, H1, h1, hy, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points,
   binned=FALSE, bgridsize, positive=FALSE, adj.positive, w, verbose=FALSE,
   nref=1e4)

\method{summary}{kroc}(object, ...)
}

\arguments{
  \item{x1,x2}{vector/matrix of data values}
  \item{H1,h1,hy}{bandwidth matrix/scalar bandwidths. If these are
    missing, \code{Hpi.kcde}, \code{hpi.kcde} is called by default.} 
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{not yet implemented}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{positive}{flag if 1-d data are positive. Default is FALSE.}
  \item{adj.positive}{adjustment applied to positive 1-d data}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{nref}{sample size used to compute reference ROC curve. Default
  is 1e3.}
  \item{object}{object of class \code{kroc}, output from \code{kroc}}
  \item{...}{other summary parameters} 
}

\value{
  A kernel ROC curve is an object of class \code{kroc} which is a list
  with fields:
  \item{x}{list of data values \code{x1, x2} - same as input}
  \item{eval.points}{points at which the estimate is evaluated}
  \item{estimate}{ROC curve estimate at \code{eval.points}}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names}
  \item{w}{weights}
  \item{tail}{"lower.tail"}
  \item{h1}{scalar bandwidth for first sample (1-d only)}
  \item{H1}{bandwidth matrix for first sample}
  \item{hy}{scalar bandwidth for ROC curve}
  \item{indices}{summary indices of ROC curve.}
}



\details{
  The usual method for computing 1-d ROC curves is based on  
  \eqn{R(z) = F_2 (\bar{F}_1^{-1}(z))}{(R(z) = F_2(bar(F)_1^(-1)(z))} where
  \eqn{\bar{F}_1, F_2}{bar(F)_1,F_2} are the respective survival and
  cumulative distribution function, and 
  \eqn{\bar{F}_1^{-1}(z)}{bar(F)_1^(-1)(z)} is the \eqn{z}{z}-th quantile of
  \eqn{\bar{F}_1}{bar(F)_1}, is not valid for multivariate data. Duong (2013),
  based on Lloyd (1998), develops an alternative formulation
  \eqn{R(z) = F_Y(z)}{R(z) = F_Y(z)} as the cumulative distribution
  function of \eqn{Y = \bar{F}_1(\bold{X}_2)}{Y=bar(F)_1(X_2)}.

  In this set-up, the values in the first sample \code{x1} should
  be larger in general that those in the second sample \code{x2}.

  If the bandwidth \code{H1} is missing from \code{kroc}, then
  the default bandwidth is the binned 2-stage plug-in selector
  \code{Hpi.kcde(, nstage=2, binned=TRUE)}. Likewise for
  missing \code{h1,hy}. A bandwidth matrix \code{H1} is required for
  \code{x1} for d>1, but the second bandwidth \code{hy} is always a 
  scalar since Y is a 1-d variable. These bandwidth selectors are optimal for cumulative
  distribution/survival functions, see Duong (2013).

  To compute the reference ROC curve, either \code{kcde} or
  \code{pmvnorm} (from the \pkg{mvtnorm} package) is used. The former is
  faster because it's a binned \code{kcde} computation, whereas the
  latter is a more exact computation of the multivariate normal CDF and
  thus can give a smoother ROC curve. 
  
  Binning/exact estimation and positive 1-d data behaviour is the same
  as for \code{\link{kde}}.
   
  --The \code{summary} method for \code{kroc} objects prints out the
  summary indices of the ROC curve, as contained in the \code{indices} field, namely
  the AUC (area under the curve) and Youden index, see Duong (2013). 

}

\references{
  Duong, T. (2013) Non-parametric kernel estimation of multivariate
  cumulative distribution functions and receiver operating
  characteristic curves. Submitted.

  Lloyd, C. (1998) Using smoothed receiver operating curves to summarize
  and compare diagnostic systems. \emph{Journal of the American
  Statistical Association}, \bold{93}, 1356--1364.
}

\seealso{\code{\link{kcde}}}
 
\examples{
samp <- 1000
x <- rnorm.mixt(n=samp, mus=0, sigmas=1, props=1)
y <- rnorm.mixt(n=samp, mus=0.5, sigmas=1, props=1)
Rhat <- kroc(x1=x, x2=y)  
summary(Rhat)
}

\keyword{ smooth }
