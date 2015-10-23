\name{Hlscv}
\alias{Hlscv}
\alias{Hlscv.diag}
\alias{Hucv}
\alias{Hucv.diag}
\alias{hlscv}
\alias{hucv}


\title{Least-squares cross-validation (LSCV) bandwidth matrix selector
  for multivariate data}
\description{
LSCV bandwidth for 1- to 6-dimensional data}
\usage{
Hlscv(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, 
      verbose=FALSE, optim.fun="nlm", trunc)
Hlscv.diag(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0, 
      verbose=FALSE, optim.fun="nlm", trunc)
hlscv(x, binned=TRUE, bgridsize, amise=FALSE, deriv.order=0)
Hucv(...)
Hucv.diag(...)
hucv(...)
}


\arguments{
  \item{x}{vector or matrix of data values}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal LSCV value. Default is FALSE.}	
  \item{deriv.order}{derivative order} 
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{optim.fun}{optimiser function: one of \code{\link{nlm}} or \code{\link{optim}}}
  \item{trunc}{parameter to control truncation for numerical
  optimisation. Default is 4 for density.deriv>0, otherwise no
  truncation. For details see below.}
  \item{...}{parameters as above}
}

\value{
LSCV bandwidth. If \code{amise=TRUE} then the minimal LSCV value is returned too. 
}

\references{
  Bowman, A. (1984) An alternative method of cross-validation for the
  smoothing of kernel density estimates. \emph{Biometrika}. \bold{71},
  353-360.

  Rudemo, M. (1982) Empirical choice of histograms and kernel density
  estimators. \emph{Scandinavian Journal of Statistics}. \bold{9},
  65-78.
}

\details{\code{hlscv} is the univariate LSCV
  selector of Bowman (1984) and Rudemo (1982). \code{Hlscv} is a
  multivariate generalisation of this. Use \code{Hlscv} for full
  bandwidth matrices and \code{Hlscv.diag} for diagonal bandwidth matrices. 
  \code{Hucv}, \code{Hucv.diag} and \code{hucv} are aliases with UCV
  unbiased cross validation instead of LSCV. 
  
  Truncation of the parameter space is usually required for the LSCV selector,
  for r > 0, to find a reasonable solution to the numerical optimisation. 
  If a candidate matrix \code{H} is
  such that \code{det(H)} is not in \code{[1/trunc, trunc]*det(H0)} or  
  \code{abs(LSCV(H)) > trunc*abs(LSCV0)} then the \code{LSCV(H)} is reset to \code{LSCV0} where 
  \code{H0=Hns(x)} and \code{LSCV0=LSCV(H0)}.

  For details about the advanced options for \code{binned,Hstart}, 
  see \code{\link{Hpi}}. 
}

\seealso{\code{\link{Hbcv}}, \code{\link{Hpi}}, \code{\link{Hscv}}}

\examples{
library(MASS)
data(forbes)
Hlscv(forbes)
hlscv(forbes$bp)
}

\keyword{ smooth }
