\name{Hpi}
\alias{Hpi}
\alias{Hpi.diag}
\alias{hpi}



\title{Plug-in bandwidth selector}

\description{
  Plug-in bandwidth for for 1- to 6-dimensional data.}

\usage{
Hpi(x, nstage=2, pilot="samse", pre="sphere", Hstart,
    binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
Hpi.diag(x, nstage=2, pilot="samse", pre="scale", Hstart,
    binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
hpi(x, nstage=2, binned=TRUE, bgridsize)
}

\arguments{
  \item{x}{vector or matrix of data values}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"amse"} = AMSE pilot bandwidths,
    \code{"samse"} = single SAMSE pilot bandwidth, 
    \code{"vamse"} = single VAMSE pilot bandwidth, 
    \code{"unconstr"} = unconstrained pilot bandwidth}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{kfold}{value for k-fold bandwidth selection. See details below.}
  \item{deriv.order}{derivative order}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
}

\value{
Plug-in bandwidth. 
If \code{amise=TRUE} then the minimal scaled PI value is returned too.}



\details{\code{hpi} is the univariate plug-in 
  selector of Wand & Jones (1994), i.e. it is exactly the same as \pkg{KernSmooth}'s \code{dpik}.
  \code{Hpi} is a multivariate generalisation of this. Use \code{Hpi} for full bandwidth matrices and \code{Hpi.diag}
  for diagonal bandwidth matrices.

  For AMSE pilot bandwidths, see Wand & Jones (1994). For
  SAMSE pilot bandwidths, see Duong & Hazelton (2003).  The latter is a
  modification of the former, in order to remove any possible problems
  with non-positive definiteness. Unconstrained pilot bandwidths are
  from Chacon & Duong (2010). VAMSE pilots are a hybrid of SAMSE and unconstrained pilots.

  For d = 1, 2, 3, 4 and \code{binned=TRUE}, 
  estimates are computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly. 
  
  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where
  \eqn{k=\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{k=4/(n*(d + 2))^(2/(d+
    4))}, n = sample size, d = dimension of data.
 
  For large samples, k-fold bandwidth selection can significantly reduce computation time. The full	
  data sample is partitioned into k sub-samples and a bandwidth matrix is computed for each of these
  sub-samples. The bandwidths are averaged and re-weighted to serve as a proxy for the full sample selector.
  (Temporarily disabled).
}

\references{
  Chacon, J.E. & Duong, T. (2010) Multivariate plug-in bandwidth
  selection with unconstrained pilot matrices. \emph{Test}, \bold{19}, 375-398.
  
  Duong, T. & Hazelton, M.L. (2003) Plug-in bandwidth matrices for
  bivariate kernel density estimation. \emph{Journal of Nonparametric
	Statistics}, \bold{15}, 17-30.
  
  Sheather, S.J. & Jones, M.C. (1991) A reliable data-based bandwidth selection
  method for kernel density estimation. \emph{Journal of the Royal
	Statistical Society, Series B}, \bold{53}, 683-690.
  
  Wand, M.P. & Jones, M.C. (1994) Multivariate plugin bandwidth
  selection. \emph{Computational Statistics}, \bold{9}, 97-116.
}
  
\examples{
data(unicef)
Hpi(unicef)
Hpi(unicef, pilot="unconstr")
Hpi.diag(unicef, binned=TRUE)
hpi(unicef[,1])
}

\keyword{ smooth }
