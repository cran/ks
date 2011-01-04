\name{Hlscv}
\alias{Hlscv}
\alias{Hlscv.diag}
\alias{hlscv}


\title{Least-squares cross-validation (LSCV) bandwidth matrix selector
  for multivariate data}
\description{
LSCV bandwidth for 1- to 6-dimensional data}
\usage{
Hlscv(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
Hlscv.diag(x, Hstart, binned=FALSE, bgridsize, amise=FALSE, kfold=1, deriv.order=0, verbose=FALSE)
hlscv(x, binned=TRUE, bgridsize, deriv.order=0)
}


\arguments{
  \item{x}{vector or matrix of data values}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal LSCV value. Default is FALSE.}	
  \item{kfold}{value for k-fold bandwidth selection. See details below.} 
  \item{deriv.order}{derivative order} 
  \item{verbose}{flag to print out progress information. Default is FALSE.}
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

\details{\code{hlscv} is the univariate SCV
  selector of Bowman (1984) and Rudemo (1982). \code{Hlscv} is a
  multivariate generalisation of this. Use \code{Hlscv} for full bandwidth matrices and \code{Hlscv.diag}
  for diagonal bandwidth matrices. 

  For details about the advanced options for \code{binned,Hstart,kfold}, 
  see \code{\link{Hpi}}. 
}

\seealso{ \code{\link{Hbcv}}, \code{\link{Hscv}}}

\examples{
library(MASS)
data(forbes)
Hlscv(forbes)
Hlscv.diag(forbes, binned=TRUE)
hlscv(forbes$bp)
}

\keyword{ smooth }
