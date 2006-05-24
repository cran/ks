\name{Hscv}
\alias{Hscv}

\title{Smoothed cross-validation (SCV) bandwidth matrix selector
  for bivariate data}
\description{
SCV bandwidth matrix for 2- to 6-dimensional data.}
\usage{
Hscv(x, pre="sphere", Hstart, binned=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{if TRUE (default) used binned estimation otherwise use
    exact kernel density estimation}
}

\value{
Full SCV bandwidth matrix.  Please note that a diagonal version of this
selector is not available. 
}


\details{ This SCV selector is a generalisation of the univariate SCV
  selector of Jones, Marron \& Park (1991).

  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.
  
  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data.}

\references{
  Jones, M.C., Marron, J.~S. \& Park, B.U. (1991) A simple root n
  bandwidth selector. \emph{The Annals of Statistics} \bold{19}, 1919--1932.

  Duong, T. \& Hazelton, M.L. (2005) \emph{Cross-validation bandwidth
    matrices for multivariate kernel density estimation}.  Scandinavian Journal 
  of Statistics. \bold{32}, 485-506.
} 

\seealso{ \code{\link{Hlscv}}, \code{\link{Hbcv}}}

\examples{
data(unicef)
Hscv(unicef)
}
\keyword{ smooth }
