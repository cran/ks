\name{plotmixt}
\alias{plotmixt}

\title{Plot for 2- and 3-dimensional normal mixture density functions.}
\description{
  Plot for 2- and 3-dimensional normal mixture density functions.
}


\usage{
plotmixt(mus, Sigmas, props, ...)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{...}{other graphics parameters}
  }

\value{
  2-d plot is sent to graphics window; 3-d plot to RGL window.
}


\details{
  See the graphics parameter options in \code{?plot.kde}.
}
 

\examples{

## bivariate example
mus <- rbind(c(0,0), c(-1,1))
Sigma <- matrix(c(1, 0.7, 0.7, 1), nr=2, nc=2) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus, Sigmas, props, ncont=6, drawlabels=FALSE)

## trivariate example
mus <- rbind(c(0,0,0), c(-1,1,1))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus, Sigmas, props, origin=c(-3,-3,-3))

}

\keyword{ hplot}
