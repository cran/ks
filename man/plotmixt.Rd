\name{plotmixt}
\alias{plotmixt}

\title{Plot for 1- to 3-dimensional normal and t-mixture density functions.}
\description{
  Plot for 1- to 3-dimensional normal and t-mixture density functions.
}


\usage{
plotmixt(mus, Sigmas, props, dfs, dist="normal", draw=TRUE, sigmas, ...)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices (2-d, 3-d)}
  \item{sigmas}{vector of standard deviations (1-d)}
  \item{props}{vector of mixing proportions}
  \item{dfs}{vector of degrees of freedom}
  \item{dist}{"normal" - normal mixture, "t" - t-mixture}
  \item{draw}{flag to draw plot. Default is TRUE.}
  \item{...}{other graphics parameters}
  }

\value{
  If \code{draw=TRUE}, the 1-d, 2-d plot is sent to graphics window, 3-d plot to RGL window.
  If \code{draw=FALSE}, then a \code{kde}-like object is returned.
}


\details{
  See the graphics parameter options in \code{?plot.kde}.
}
 

\examples{
## unvariate 
mus <- c(0,1)
sigmas <- c(0.2, 0.5) 
props <- c(1/2, 1/2)
plotmixt(mus=mus, sigmas=sigmas, props=props)

## bivariate 
mus <- rbind(c(0,0), c(-1,1))
Sigma <- matrix(c(1, 0.7, 0.7, 1), nr=2, nc=2) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus, Sigmas, props)

## trivariate 
mus <- rbind(c(0,0,0), c(-1,0.5,1.5))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus, Sigmas, props, dfs=c(3,8), dist="t")
}

\keyword{ hplot}
