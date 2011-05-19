\name{Hscv}
\alias{Hscv}
\alias{Hscv.diag}
\alias{hscv}
\title{Smoothed cross-validation (SCV) bandwidth selector}

\description{
  SCV bandwidth for 1- to 6-dimensional data.}

\usage{
Hscv(x, nstage=2, pre="sphere", pilot="samse", Hstart, binned=FALSE, 
     bgridsize, amise=FALSE, kfold=1, verbose=FALSE, optim.fun="nlm")
Hscv.diag(x, nstage=2, pre="scale", pilot="samse", Hstart, binned=FALSE, 
     bgridsize, amise=FALSE, kfold=1, verbose=FALSE, optim.fun="nlm")
hscv(x, nstage=2, binned=TRUE, bgridsize, plot=FALSE)
}
\arguments{
  \item{x}{vector or matrix of data values}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} = pre-sphering}
  \item{pilot}{\code{"amse"} = AMSE pilot bandwidths,
    \code{"samse"} = single SAMSE pilot bandwidth,
    \code{"vamse"} = single VAMSE pilot bandwidth, 
    \code{"unconstr"} = unconstrained pilot bandwidth matrix}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal scaled SCV value. Default is FALSE.}
  \item{kfold}{value for k-fold bandwidth selection. See details below.} 
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{optim.fun}{optimiser function: one of \code{\link{nlm}} or \code{\link{optim}}.}
  \item{nstage}{number of stages in the SCV bandwidth selector (1 or 2)}
  \item{plot}{flag to display plot of SCV(h) vs h (1-d only). Default is FALSE.}
}

\value{
SCV bandwidth. If \code{amise=TRUE} then the minimal scaled SCV value is returned too. 
}


\details{ \code{hsv} is the univariate SCV
  selector of Jones, Marron & Park (1991). \code{Hscv} is a
  multivariate generalisation of this, see Duong \& Hazelton (2005).
  Use \code{Hscv} for full bandwidth matrices and \code{Hscv.diag}
  for diagonal bandwidth matrices. 
 
  For AMSE pilot bandwidths, see Wand & Jones (1994). For
  SAMSE pilot bandwidths, see Duong & Hazelton (2003).  The latter is a
  modification of the former, in order to remove any possible problems
  with non-positive definiteness. Unconstrained pilot bandwidths are
  from Chacon & Duong (2011). VAMSE pilots are a hybrid of SAMSE and unconstrained pilots.

  For d = 1, the selector \code{hscv} is not always stable for large
  sample sizes with binning.
  Examine the plot from \code{hscv(, plot=TRUE)} to
  determine the appropriate smoothness of the SCV function. Any
  non-smoothness is due to the discretised nature of binned estimation.
  
  For details about the advanced options for \code{amise, binned, Hstart, kfold}, 
  see \code{\link{Hpi}}.
}

\references{
  Chacon, J.E. & Duong, T. (2011) Unconstrained pilot selectors for smoothed cross
   validation. Under revision.
  
  Duong, T. & Hazelton, M.L. (2003) Plug-in bandwidth matrices for
  bivariate kernel density estimation. \emph{Journal of Nonparametric
	Statistics}, \bold{15}, 17-30.
  
  Duong, T. & Hazelton, M.L. (2005) Cross-validation bandwidth
    matrices for multivariate kernel density estimation. \emph{Scandinavian Journal 
  of Statistics}. \bold{32}, 485-506.

  Wand, M.P. & Jones, M.C. (1994) Multivariate plugin bandwidth
  selection. \emph{Computational Statistics}, \bold{9}, 97-116. 
  Jones, M.C., Marron, J.S. & Park, B.U. (1991) A simple root n
  bandwidth selector. \emph{Annals of Statistics} \bold{19}, 1919-1932.

} 

\seealso{ \code{\link{Hlscv}}, \code{\link{Hbcv}}, \code{\link{Hpi}}}

\examples{
data(unicef)
Hscv(unicef)
Hscv.diag(unicef, binned=TRUE)
hscv(unicef[,1])
}

\keyword{ smooth }
