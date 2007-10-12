\name{Hpi, Hpi.diag}
\alias{Hpi}
\alias{Hpi.diag}

\title{Plug-in bandwidth matrix selector for multivariate data}
\description{
Plug-in bandwidth matrix for for 2- to 6-dimensional data.}
\usage{
Hpi(x, nstage=2, pilot="samse", pre="sphere", Hstart,
    binned=FALSE, bgridsize)
Hpi.diag(x, nstage=2, pilot="amse", pre="scale", Hstart,
    binned=FALSE, bgridsize)
}

\arguments{
  \item{x}{matrix of data values}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"amse"} = AMSE pilot bandwidths,
    \code{"samse"} = single SAMSE pilot bandwidth}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - required only if binned=TRUE}
}

\value{
Plug-in bandwidth matrix.
}
\references{ Wand, M.P. \& Jones, M.C. (1994) Multivariate plugin bandwidth
    selection.  \emph{Computational Statistics}, \bold{9}, 97-116.
  
  Duong, T. \& Hazelton, M.L. (2003) Plug-in bandwidth matrices for
    bivariate kernel density estimation. \emph{Journal of Nonparametric
  Statistics}, \bold{15}, 17-30.
}

\details{ Use \code{Hpi} for full bandwidth matrices and \code{Hpi.diag}
  for diagonal bandwidth matrices.

  For AMSE pilot bandwidths, see Wand \& Jones (1994). For
  SAMSE pilot bandwidths, see Duong \& Hazelton (2003).  The latter is a
  modification of the former, in order to remove any possible problems
  with non-positive definiteness. 

  For d = 1, 2, 3, 4 and \code{binned=TRUE}, 
  the density estimate is computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.
  
  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k =
  \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+
    4))}, n = sample size, d = dimension of data.
}


\examples{
## bivariate example
data(unicef)
Hpi(unicef, nstage=2, pilot="samse", pre="sphere")
Hpi(unicef, nstage=2, pilot="samse", pre="sphere", binned=TRUE)

\dontrun{
## 4-variate example
x <- rmvnorm.mixt(n=1000, mus=rep(0,4), diag(4))
Hpi.diag(x, nstage=1, pilot="samse", pre="scale") 
Hpi.diag(x, nstage=1, pilot="samse", pre="scale", binned=TRUE)
} 
}

\keyword{ smooth }
