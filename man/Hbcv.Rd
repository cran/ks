\name{Hbcv, Hbcv.diag}
\alias{Hbcv}
\alias{Hbcv.diag}

\title{Biased cross-validation (BCV) bandwidth matrix selector
  for bivariate data}
\description{
BCV bandwidth matrix for bivariate data.}
\usage{
Hbcv(x, whichbcv=1, Hstart)
Hbcv.diag(x, whichbcv=1, Hstart)
}

\arguments{
  \item{x}{matrix of data values}
  \item{whichbcv}{1 = BCV1, 2 = BCV2.  See details below}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}  
}

\value{
BCV bandwidth matrix.   
}
\references{Sain, S.R, Baggerly, K.A. \& Scott, D.W. (1994)
  Cross-validation of multivariate densities. \emph{Journal of the
  American Statistical Association}. \bold{82}, 1131-1146.
  
  Duong, T. \& Hazelton, M.L. (2005) Cross-validation bandwidth
    matrices for multivariate kernel density estimation. \emph{Scandinavian Journal 
    of Statistics}. \bold{32}, 485-506.
  
}

\details{
  Use \code{Hbcv} for full bandwidth matrices and \code{Hbcv.diag}
  for diagonal bandwidth matrices. These selectors are only
  available for bivariate data.

  There are two types of BCV criteria
  considered here.  They are known as BCV1 and BCV2, from Sain, Baggerly
  \& Scott (1994) and they only differ slightly. These BCV
  surfaces can have multiple minima and so it can be quite difficult to
  locate the most appropriate minimum.

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data.
}

\seealso{ \code{\link{Hlscv}}, \code{\link{Hscv}}}

\note{
  It can be difficult to find an appropriate (local) minimum of the BCV
  criterion. Some times, there can be no local minimum at all so there
  may be no finite BCV selector.
}

\examples{
data(unicef)
Hbcv(unicef)
Hbcv.diag(unicef)
}
\keyword{ smooth }
