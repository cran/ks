\name{Hlscv}
\alias{Hlscv}
\alias{Hlscv.diag}


\title{Least-squares cross-validation (LSCV) bandwidth matrix selector
  for bivariate data}
\description{
LSCV bandwidth matrix for bivariate data.}
\usage{
Hlscv(x, Hstart)
Hlscv.diag(x, Hstart)
}

\arguments{
  \item{x}{matrix of data values}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation} 
}

\value{
LSCV bandwidth matrix.
}
\references{Sain, S.R, Baggerly, K.A \& Scott, D.W. (1994)
  \emph{Cross-validation of multivariate densities}. Journal of the
  American Statistical Association. \bold{82}, 1131-1146.   

  Duong, T. \& Hazelton, M.L. (2004) \emph{Cross-validation bandwidth
    matrices for multivariate kernel density estimation}.  Scandinavian Journal 
    of Statistics. In press. 
}

\details{Use \code{Hlscv} for full bandwidth matrices and \code{Hlscv.diag}
  for diagonal bandwidth matrices.

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data. }

\seealso{ \code{\link{Hbcv}}, \code{\link{Hscv}}}

\examples{
data(unicef)
Hlscv(unicef)
Hlscv.diag(unicef)
}
\keyword{ smooth }
