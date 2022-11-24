\name{kda}
\alias{Hkda}
\alias{Hkda.diag}
\alias{kda}
\alias{hkda}
\alias{predict.kda}
\alias{compare}
\alias{compare.kda.diag.cv}
\alias{compare.kda.cv}
\title{Kernel discriminant analysis (kernel classification)}

\description{
  Kernel discriminant analysis (kernel classification) for 1- to d-dimensional data.
}

\usage{
kda(x, x.group, Hs, hs, prior.prob=NULL, gridsize, xmin, xmax, supp=3.7,
  eval.points, binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE,
  kde.flag=TRUE)
Hkda(x, x.group, Hstart, bw="plugin", ...)
Hkda.diag(x, x.group, bw="plugin", ...)
hkda(x, x.group, bw="plugin", ...)

\method{predict}{kda}(object, ..., x)

compare(x.group, est.group, by.group=FALSE)
compare.kda.cv(x, x.group, bw="plugin", prior.prob=NULL, Hstart, by.group=FALSE,
   verbose=FALSE, recompute=FALSE, ...)
compare.kda.diag.cv(x, x.group, bw="plugin", prior.prob=NULL, by.group=FALSE, 
   verbose=FALSE, recompute=FALSE, ...)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{Hs,hs}{(stacked) matrix of bandwidth matrices/vector of scalar
  bandwidths. If these are missing, \code{Hkda} or \code{hkda} is called by default.}
  \item{prior.prob}{vector of prior probabilities}
  \item{gridsize}{vector of grid sizes}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{vector or matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights. Not yet implemented.}
  \item{compute.cont}{flag for computing 1\% to 99\% probability contour levels. Default is TRUE.}
  \item{approx.cont}{flag for computing approximate probability contour
    levels. Default is TRUE.}
  \item{kde.flag}{flag for computing KDE on grid. Default is TRUE.}
  \item{object}{object of class \code{kda}}
  \item{bw}{bandwidth: "plugin" = plug-in, "lscv" = LSCV, 
    "scv" = SCV}
  \item{Hstart}{(stacked) matrix of initial bandwidth matrices, used in
    numerical optimisation}
  \item{est.group}{vector of estimated group labels}
  \item{by.group}{flag to give results also within each group}
  \item{verbose}{flag for printing progress information. Default is FALSE.}
  \item{recompute}{flag for recomputing the bandwidth matrix after
    excluding the i-th data item}
  \item{...}{other optional parameters for bandwidth selection, see
    \code{\link{Hpi}}, \code{\link{Hlscv}}, \code{\link{Hscv}}} 
}

\value{
  --For \code{kde.flag=TRUE}, a kernel discriminant analysis is an object of class \code{kda} which is a list with fields
  \item{x}{list of data points, one for each group label}
  \item{estimate}{list of density estimates at \code{eval.points}, one for each group label}  
  \item{eval.points}{vector or list of points that the estimate is evaluated at, one
    for each group label}
  \item{h}{vector of bandwidths (1-d only)}
  \item{H}{stacked matrix of bandwidth matrices or vector of bandwidths}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{w}{vector of weights}
  \item{prior.prob}{vector of prior probabilities}
  \item{x.group}{vector of group labels - same as input}
  \item{x.group.estimate}{vector of estimated group labels. If the test data
    \code{eval.points} are given then these are classified. Otherwise
    the training data \code{x} are classified.}
  For \code{kde.flag=FALSE}, which is always the case for \eqn{d>3}{d>3},
  then only the vector of estimated group labels is returned. 
  
  --The result from \code{Hkda} and \code{Hkda.diag} is a stacked matrix
  of bandwidth matrices, one for each training data group. The result
  from \code{hkda} is a vector of bandwidths, one for each training group.  

  --The \code{compare} functions create a comparison between the true
  group labels \code{x.group} and the estimated ones. 
  It returns a list with fields
  \item{cross}{cross-classification table with the rows
    indicating the true group and the columns the estimated group}
  \item{error}{misclassification rate (MR)}
    
  In the case where the test data are independent of the
  training data, \code{compare} computes MR = (number of points wrongly
  classified)/(total number of points). In the case where the test data
  are not independent e.g. 
  we are classifying the training data set itself, then the cross
  validated estimate of MR is more appropriate. These
  are implemented as \code{compare.kda.cv} (unconstrained bandwidth
  selectors) and \code{compare.kda.diag.cv} (for diagonal bandwidth
  selectors). These functions are only available for d > 1.

  If \code{by.group=FALSE} then only the total MR rate is given. If it
  is set to TRUE, then the MR rates for each class are also given
  (estimated number in group divided by true number).
}

\references{
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York
}

\details{
  If the bandwidths \code{Hs} are missing from \code{kda}, then the
  default bandwidths are the plug-in selectors \code{Hkda(, bw="plugin")}.
  Likewise for missing \code{hs}. Valid options for \code{bw}
  are \code{"plugin"}, \code{"lscv"} and \code{"scv"} which in turn call
  \code{\link{Hpi}}, \code{\link{Hlscv}} and \code{\link{Hscv}}.

  The effective support, binning, grid size, grid range, positive
  parameters are the same as \code{\link{kde}}.
  
  If prior probabilities are known then set \code{prior.prob} to these.
  Otherwise \code{prior.prob=NULL} uses the sample
  proportions as estimates of the prior probabilities.
  
  For \pkg{ks} \eqn{\geq}{>=} 1.8.11, \code{kda.kde} has been subsumed
  into \code{kda}, so all prior calls to \code{kda.kde} can be replaced
  by \code{kda}. To reproduce the previous behaviour of \code{kda}, the
  command is \code{kda(, kde.flag=FALSE)}.
}

\seealso{\code{\link{plot.kda}}}

\examples{
set.seed(8192)
x <- c(rnorm.mixt(n=100, mus=1), rnorm.mixt(n=100, mus=-1))
x.gr <- rep(c(1,2), times=c(100,100))
y <- c(rnorm.mixt(n=100, mus=1), rnorm.mixt(n=100, mus=-1))
y.gr <- rep(c(1,2), times=c(100,100))
kda.gr <- kda(x, x.gr)
y.gr.est <- predict(kda.gr, x=y)
compare(y.gr, y.gr.est)

## See other examples in ? plot.kda
}

\keyword{smooth}
