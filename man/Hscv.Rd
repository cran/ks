\name{Hscv}
\alias{Hscv}
\alias{Hscv.diag}
\alias{hscv}
\title{Smoothed cross-validation (SCV) bandwidth selector}

\description{
  SCV bandwidth for 1- to 6-dimensional data.}

\usage{
Hscv(x, nstage=2, pre="sphere", pilot, Hstart, binned, 
     bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
Hscv.diag(x, nstage=2, pre="scale", pilot, Hstart, binned, 
     bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
hscv(x, nstage=2, binned=TRUE, bgridsize, plot=FALSE)
}

\arguments{
  \item{x}{vector or matrix of data values}
  \item{pre}{"scale" = \code{\link{pre.scale}}, "sphere" = \code{\link{pre.sphere}}}
  \item{pilot}{"amse" = AMSE pilot bandwidths \cr
    "samse" = single SAMSE pilot bandwidth \cr 
    "unconstr" = single unconstrained pilot bandwidth \cr
    "dscalar" = single pilot bandwidth for deriv.order>0 \cr
    "dunconstr" = single unconstrained pilot bandwidth for deriv.order>0}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal scaled SCV value. Default is FALSE.}
  \item{deriv.order}{derivative order}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{optim.fun}{optimiser function: one of \code{nlm} or \code{optim}}
  \item{nstage}{number of stages in the SCV bandwidth selector (1 or 2)}
  \item{plot}{flag to display plot of SCV(h) vs h (1-d only). Default is FALSE.}
}

\value{
SCV bandwidth. If \code{amise=TRUE} then the minimal scaled SCV value is returned too. 
}

\details{ \code{hscv} is the univariate SCV
  selector of Jones, Marron & Park (1991). \code{Hscv} is a
  multivariate generalisation of this, see Duong & Hazelton (2005).
  Use \code{Hscv} for unconstrained bandwidth matrices and \code{Hscv.diag}
  for diagonal bandwidth matrices. 
 
  The default pilot is \code{"samse"} for d=2, r=0, and
  \code{"dscalar"} otherwise. For SAMSE pilot bandwidths, see Duong &
  Hazelton (2005). Unconstrained and higher order derivative pilot
  bandwidths are from Chacon & Duong (2011).  

  For d=1, the selector \code{hscv} is not always stable for large
  sample sizes with binning.
  Examine the plot from \code{hscv(, plot=TRUE)} to
  determine the appropriate smoothness of the SCV function. Any
  non-smoothness is due to the discretised nature of binned estimation.

  For details about the advanced options for \code{binned, Hstart, optim.fun}, 
  see \code{\link{Hpi}}.
}

\references{
  Chacon, J.E. & Duong, T. (2011) Unconstrained pilot selectors for
  smoothed cross validation. \emph{Australian & New Zealand Journal of
  Statistics}, \bold{53}, 331-351. 

  Duong, T. & Hazelton, M.L. (2005) Cross-validation bandwidth
  matrices for multivariate kernel density
  estimation. \emph{Scandinavian Journal of Statistics}, \bold{32}, 485-506.
   
  Jones, M.C., Marron, J.S. & Park, B.U. (1991) A simple root \eqn{n}{n}
  bandwidth selector. \emph{Annals of Statistics}, \bold{19}, 1919-1932.
} 

\seealso{\code{\link{Hbcv}}, \code{\link{Hlscv}}, \code{\link{Hpi}}}

\examples{
data(unicef)
Hscv(unicef)
hscv(unicef[,1])
}

\keyword{smooth}
