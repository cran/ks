\name{Hpi}
\alias{Hpi}
\alias{Hpi.diag}
\alias{hpi}
\title{Plug-in bandwidth selector}

\description{
  Plug-in bandwidth for for 1- to 6-dimensional data.}

\usage{
Hpi(x, nstage=2, pilot, pre="sphere", Hstart, binned, bgridsize,
   amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
Hpi.diag(x, nstage=2, pilot, pre="scale", Hstart, binned, bgridsize,
   amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
hpi(x, nstage=2, binned=TRUE, bgridsize, deriv.order=0)
}

\arguments{
  \item{x}{vector or matrix of data values}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{"amse" = AMSE pilot bandwidths \cr
    "samse" = single SAMSE pilot bandwidth \cr 
    "unconstr" = single unconstrained pilot bandwidth \cr
    "dscalar" = single pilot bandwidth for deriv.order >= 0 \cr
    "dunconstr" = single unconstrained pilot bandwidth for deriv.order >= 0}
  \item{pre}{"scale" = \code{\link{pre.scale}}, "sphere" = \code{\link{pre.sphere}}}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{deriv.order}{derivative order}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{optim.fun}{optimiser function: one of \code{nlm} or \code{optim}}
}

\value{
Plug-in bandwidth. 
If \code{amise=TRUE} then the minimal scaled PI value is returned too.}

\details{\code{hpi(,deriv.order=0)} is the univariate plug-in 
  selector of Wand & Jones (1994), i.e. it is exactly the same as
  \pkg{KernSmooth}'s \code{dpik}. For deriv.order>0, the formula is
  taken from Wand & Jones (1995). \code{Hpi} is a multivariate
  generalisation of this. Use \code{Hpi} for unconstrained bandwidth matrices and
  \code{Hpi.diag} for diagonal bandwidth matrices.

  The default pilot is \code{"samse"} for d=2,r=0, and
  \code{"dscalar"} otherwise.
  For AMSE pilot bandwidths, see Wand & Jones (1994). For
  SAMSE pilot bandwidths, see Duong & Hazelton (2003).  The latter is a
  modification of the former, in order to remove any possible problems
  with non-positive definiteness. Unconstrained and higher order
  derivative pilot bandwidths are from Chacon & Duong (2010). 

  For d=1, 2, 3, 4 and \code{binned=TRUE}, 
  estimates are computed over a binning grid defined 
  by \code{bgridsize}. Otherwise it's computed exactly.  
  If \code{Hstart} is not given then it defaults to \code{Hns(x)}.

  For \pkg{ks} \eqn{\geq}{>=} 1.11.1, the default optimisation function is
  \code{optim.fun="optim"}. To reinstate the previous functionality, use
  \code{optim.fun="nlm"}.  
}

\references{
  Chacon, J.E. & Duong, T. (2010) Multivariate plug-in bandwidth
  selection with unconstrained pilot matrices. \emph{Test}, \bold{19}, 375-398.
  
  Duong, T. & Hazelton, M.L. (2003) Plug-in bandwidth matrices for
  bivariate kernel density estimation. \emph{Journal of Nonparametric
    Statistics}, \bold{15}, 17-30.
  
  Sheather, S.J. & Jones, M.C. (1991) A reliable data-based bandwidth selection
  method for kernel density estimation. \emph{Journal of the Royal
    Statistical Society Series B}, \bold{53}, 683-690.
  
  Wand, M.P. & Jones, M.C. (1994) Multivariate plug-in bandwidth
  selection. \emph{Computational Statistics}, \bold{9}, 97-116.
}

\seealso{\code{\link{Hbcv}}, \code{\link{Hlscv}}, \code{\link{Hscv}}}
  
\examples{
data(unicef)
Hpi(unicef, pilot="dscalar")
hpi(unicef[,1])
}

\keyword{smooth}
