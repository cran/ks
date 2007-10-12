\name{Hscv,Hscv.diag}
\alias{Hscv}
\alias{Hscv.diag}

\title{Smoothed cross-validation (SCV) bandwidth matrix selector
  for 2- to 6-dimensional data}
\description{
SCV bandwidth matrix for 2- to 6-dimensional data.}
\usage{
Hscv(x, pre="sphere", Hstart, binned=FALSE, bgridsize)
Hscv.diag(x, pre="scale", Hstart, binned=FALSE, bgridsize)
}

\arguments{
  \item{x}{matrix of data values}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - 
  required only if \code{binned=TRUE}}
}

\value{
SCV bandwidth matrix.  
}


\details{ This SCV selector is a generalisation of the univariate SCV
  selector of Jones, Marron \& Park (1991).

  For d = 1, 2, 3, 4 and \code{binned=TRUE}, 
  the density estimate is computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.
  
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

\seealso{ \code{\link{Hlscv}}, \code{\link{Hbcv}}}

\examples{
data(unicef)
Hscv(unicef)
Hscv(unicef, binned=TRUE)
Hscv.diag(unicef)
}
\keyword{ smooth }
