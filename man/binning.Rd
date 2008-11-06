\name{binning}
\alias{binning}

\title{Linear binning for multivariate data}
\description{
  Linear binning for 1- to 4-dimensional data.
}
\usage{
binning(x, H, h, bgridsize, xmin, xmax, supp=3.7)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{h}{scalar bandwidth}
  \item{xmin}{vector of minimum values for grid}
  \item{xmax}{vector of maximum values for grid}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{bgridsize}{vector of binning grid sizes}
}
  

\value{
  Returns a list with 2 fields
  \item{counts}{linear binning counts}
  \item{eval.points}{vector (d=1) or list (d>2) of grid points in each
	dimension }
}

\details{
  Code is used courtesy of Matt Wand. Default \code{bgridsize} are:
  401 (d=1), rep(151, 2), rep(51, 3), rep(21,4). 
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\examples{
data(unicef)
binning(x=unicef)
binning(x=unicef, xmin=c(0, 20), xmax=c(350, 100), H=diag(2), bgridsize=c(201,201))
}

\keyword{smooth}