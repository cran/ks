\name{binning}
\alias{binning}

\title{Linear binning for multivariate data}
\description{
  Linear binning for 1- to 4-dimensional data.
}
\usage{
binning(x, H, h, bgridsize, xmin, xmax, supp=3.7, w)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix, scalar bandwidth}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal is [\code{-supp,supp}]}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights (non-negative and sum is equal to sample size). Default is a vector of all ones.}
}
  

\value{
  Returns a list with 2 fields
  \item{counts}{linear binning counts}
  \item{eval.points}{vector (d=1) or list (d>2) of grid points in each
	dimension }
}

\details{
  Code is used courtesy of M.P. Wand. Default \code{bgridsize} are
  d=1: 401; d=2: rep(151, 2); d=3: rep(31, 3); d=4: rep(21,4). 
}

\references{ Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman & Hall. London. } 

\examples{
data(unicef)
ubinned <- binning(x=unicef)
ubinned <- binning(x=unicef, xmin=c(0, 20), xmax=c(350, 100))
}

\keyword{smooth}
