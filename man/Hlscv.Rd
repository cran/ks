\name{Hlscv, Hlscv.diag}
\alias{Hlscv}
\alias{Hlscv.diag}


\title{Least-squares cross-validation (LSCV) bandwidth matrix selector
  for multivariate data}
\description{
LSCV bandwidth matrix for 2- to 6-dimensional data}
\usage{
Hlscv(x, Hstart)
Hlscv.diag(x, Hstart, binned=FALSE, bgridsize)
}

\arguments{
  \item{x}{matrix of data values}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - 
  required only if \code{binned=TRUE}}
}

\value{
LSCV bandwidth matrix.
}
\references{Sain, S.R, Baggerly, K.A \& Scott, D.W. (1994)
  Cross-validation of multivariate densities. \emph{Journal of the
  American Statistical Association}. \bold{82}, 1131-1146.   

  Duong, T. \& Hazelton, M.L. (2005) Cross-validation bandwidth
    matrices for multivariate kernel density estimation.  \emph{Scandinavian Journal 
    of Statistics}. \bold{32}, 485-506.
}

\details{Use \code{Hlscv} for full bandwidth matrices and \code{Hlscv.diag}
  for diagonal bandwidth matrices.
  
  For d = 1, 2, 3, 4 and \code{binned=TRUE}, 
  the density estimate is computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.
  
  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data. }

\seealso{ \code{\link{Hbcv}}, \code{\link{Hscv}}}

\examples{
mus <- rbind(c(-1/2,0), c(1/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/8, 1/16), c(1/16, 1/8)))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(1000, mus, Sigmas, props)
Hlscv(x)
Hlscv.diag(x)
Hlscv.diag(x, binned=TRUE)
}
\keyword{ smooth }
