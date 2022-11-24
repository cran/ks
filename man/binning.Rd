\name{binning}
\alias{binning}
\title{Linear binning for multivariate data}

\description{
  Linear binning for 1- to 4-dimensional data.
}

\usage{
binning(x, H, h, bgridsize, xmin, xmax, supp=3.7, w, gridtype="linear")
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix, scalar bandwidth}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal is [-supp,supp]}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{gridtype}{not yet implemented}
}
  
\value{
  Returns a list with 2 fields
  \item{counts}{linear binning counts}
  \item{eval.points}{vector (d=1) or list (d>=2) of grid points in each
	dimension }
}

\details{
  For \pkg{ks} \eqn{\geq}{>=} 1.10.0, binning is available for unconstrained
  (non-diagonal) bandwidth matrices. Code is used courtesy of A. &
  J. Gramacki, and M.P. Wand. Default \code{bgridsize} are 
  d=1: 401; d=2: rep(151, 2); d=3: rep(51, 3); d=4: rep(21, 4). 
}

\references{
  Gramacki, A. & Gramacki, J. (2016) FFT-based fast computation of
  multivariate kernel estimators with unconstrained bandwidth
  matrices. \emph{Journal of Computational & Graphical Statistics},
  \bold{26}, 459-462.

  Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman & Hall. London.
} 

\examples{
data(unicef)
ubinned <- binning(x=unicef)
}

\keyword{algebra}
