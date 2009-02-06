\name{Hscv, Hscv.diag, hscv}
\alias{Hscv}
\alias{Hscv.diag}
%\alias{hscv}
\title{Smoothed cross-validation (SCV) bandwidth selector}

\description{
  SCV bandwidth for 1- to 6-dimensional data.}

\usage{
Hscv(x, pre="sphere", Hstart, binned=FALSE, bgridsize)
Hscv.diag(x, pre="scale", Hstart, binned=FALSE, bgridsize)
}
%hscv(x, nstage=2, binned=TRUE, bgridsize)

\arguments{
  \item{x}{vector or matrix of data values}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - 
	required only if \code{binned=TRUE}}
  %\item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
}

\value{
SCV bandwidth.  
}


\details{ \code{hsv} is the univariate SCV
  selector of Jones, Marron \& Park (1991). \code{Hscv} is a
  multivariate generalisation of this.

  %For d = 1, the selector \code{hscv} is always computed as binned
  %estimator.
  For d = 2, 3, 4 and \code{binned=TRUE}, the density
  estimate is computed over a binning grid defined by
  \code{bgridsize}. Otherwise it's computed exactly. 
  
  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.
  
  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data.}

\references{
  Jones, M.C., Marron, J.~S. \& Park, B.U. (1991) A simple root n
  bandwidth selector. \emph{Annals of Statistics} \bold{19}, 1919--1932.

  Duong, T. \& Hazelton, M.L. (2005) Cross-validation bandwidth
    matrices for multivariate kernel density estimation. \emph{Scandinavian Journal 
  of Statistics}. \bold{32}, 485-506.
} 

\seealso{ \code{\link{Hlscv}}, \code{\link{Hbcv}}, \code{\link{Hpi}}}

\examples{
data(unicef)
Hscv(unicef)
Hscv.diag(unicef, binned=TRUE)
}
%hscv(unicef[,1])
\keyword{ smooth }
