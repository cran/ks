\name{Hlscv, Hlscv.diag, hlscv}
\alias{Hlscv}
\alias{Hlscv.diag}
\alias{hlscv}


\title{Least-squares cross-validation (LSCV) bandwidth matrix selector
  for multivariate data}
\description{
LSCV bandwidth for 1- to 6-dimensional data}
\usage{
Hlscv(x, Hstart)
Hlscv.diag(x, Hstart, binned=FALSE, bgridsize)
hlscv(x, binned=TRUE, bgridsize)
}


\arguments{
  \item{x}{vector or matrix of data values}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - 
  required only if \code{binned=TRUE}}
}

\value{
LSCV bandwidth.
}
\references{
  Bowman, A. (1984) An alternative method of cross-validation for the
  smoothing of kernel density estimates. \emph{Biometrika}. \bold{71},
  353-360.
  
  Duong, T. \& Hazelton, M.L. (2005) Cross-validation
  bandwidth matrices for multivariate kernel density estimation.
  \emph{Scandinavian Journal of Statistics}. \bold{32}, 485-506.

  Rudemo, M. (1982) Empirical choice of histograms and kernel density
  estimators. \emph{Scandinavian Journal of Statistics}. \bold{9},
  65-78.
  
  Sain, S.R, Baggerly, K.A \& Scott, D.W. (1994)
  Cross-validation of multivariate densities. \emph{Journal of the
  American Statistical Association}. \bold{82}, 1131-1146.   
}

\details{\code{hlscv} is the univariate SCV
  selector of Bowman (1984) and Rudemo (1982). \code{Hlscv} is a
  multivariate generalisation of this.

  Use \code{Hlscv} for full bandwidth matrices and \code{Hlscv.diag}
  for diagonal bandwidth matrices.
  
  For d = 1, 2, 3, 4 and \code{binned=TRUE}, 
  estimates are computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.
  
  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data. }

\seealso{ \code{\link{Hbcv}}, \code{\link{Hscv}}}

\examples{
library(MASS)
data(forbes)
Hlscv(forbes)
Hlscv.diag(forbes, binned=TRUE)
hlscv(forbes$bp)
}

\keyword{ smooth }
