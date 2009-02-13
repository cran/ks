\name{Hpi, Hpi.diag, hpi}
\alias{Hpi}
\alias{Hpi.diag}
\alias{hpi}


\title{Plug-in bandwidth selector}

\description{
  Plug-in bandwidth for for 1- to 6-dimensional data.}

\usage{
Hpi(x, nstage=2, pilot="samse", pre="sphere", Hstart,
    binned=FALSE, bgridsize, amise=FALSE)
Hpi.diag(x, nstage=2, pilot="amse", pre="scale", Hstart,
    binned=FALSE, bgridsize)
hpi(x, nstage=2, binned=TRUE, bgridsize)
}

\arguments{
  \item{x}{vector or matrix of data values}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"amse"} = AMSE pilot bandwidths,
    \code{"samse"} = single SAMSE pilot bandwidth,
	\code{"unconstr"} = unconstrained pilot bandwidth matrix}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - required only if binned=TRUE}
  \item{amise}{flag for returning estimated AMISE}
  %\item{...}{other parameters as for \code{dpik} from \pkg{KernSmooth}}
}

\value{
Plug-in bandwidth. If \code{amise=TRUE} then the plug-in
bandwidth plus the estimated AMISE is returned in a list.}



\details{ \code{hpi} is the univariate plug-in 
  selector of Sheather \& Jones (1991). \code{Hpi} is a
  multivariate generalisation of this.

  Use \code{Hpi} for full bandwidth matrices and \code{Hpi.diag}
  for diagonal bandwidth matrices.
  For AMSE pilot bandwidths, see Wand \& Jones (1994). For
  SAMSE pilot bandwidths, see Duong \& Hazelton (2003).  The latter is a
  modification of the former, in order to remove any possible problems
  with non-positive definiteness. Unconstrained pilot bandwidths are
  available for d = 1, ..., 5 (but are extremely computationally
  intensive for the latter dimensions). See Chac\'on \& Duong (2008).

  For d = 1, the selector \code{hpi} is exactly the same as
  \pkg{KernSmooth}'s \code{dpik}. This is always computed as binned
  estimator. For d = 2, 3, 4 and \code{binned=TRUE}, 
  estimates are computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.
  
  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k =
  \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+
    4))}, n = sample size, d = dimension of data.
}

\references{
  Chac\'on, J.E. \& Duong, T. (2008) Multivariate plug-in bandwidth
  selection with unconstrained pilot matrices. \emph{Submitted.}
  
  Duong, T. \& Hazelton, M.L. (2003) Plug-in bandwidth matrices for
  bivariate kernel density estimation. \emph{Journal of Nonparametric
	Statistics}, \bold{15}, 17-30.
  
  Sheather, S.J. \& Jones, M.C. (1991) A reliable data-based bandwidth selection
  method for kernel density estimatio. \emph{Journal of the Royal
	Statistical Society, Series B}, \bold{53}, 683-690.
  
  Wand, M.P. \& Jones, M.C. (1994) Multivariate plugin bandwidth
  selection.  \emph{Computational Statistics}, \bold{9}, 97-116.
}
  
\examples{
data(unicef)
Hpi(unicef)
Hpi(unicef, pilot="unconstr")
Hpi.diag(unicef, binned=TRUE)
hpi(unicef[,1])

x <- rmvnorm.mixt(n=100, rep(0,4), diag(4)); Hpi(x)
}

\keyword{ smooth }
