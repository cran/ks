\name{plotmixt}
\alias{plotmixt}
\title{Plot for 1- to 3-dimensional normal and t-mixture density functions}

\description{
  Plot for 1- to 3-dimensional normal and t-mixture density functions.
}

\usage{
plotmixt(mus, sigmas, Sigmas, props, dfs, dist="normal", draw=TRUE,
   deriv.order=0, which.deriv.ind=1, binned=TRUE, ...)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{sigmas}{vector of standard deviations (1-d)}
  \item{Sigmas}{(stacked) matrix of variance matrices (2-d, 3-d)}
  \item{props}{vector of mixing proportions}
  \item{dfs}{vector of degrees of freedom}
  \item{dist}{"normal" - normal mixture, "t" - t-mixture}
  \item{draw}{flag to draw plot. Default is TRUE.}
  \item{deriv.order}{derivative order}
  \item{which.deriv.ind}{index of which partial derivative to plot}
  \item{binned}{flag for binned estimation of contour levels. Default is TRUE.} 
  \item{...}{other graphics parameters, see \code{\link{plot.kde}}}
  }

\value{
  If \code{draw=TRUE}, the 1-d, 2-d plot is sent to graphics window, 3-d plot to 
  graphics/RGL window. If \code{draw=FALSE}, then a \code{kdde}-like object is returned.
}

\examples{
## bivariate 
mus <- rbind(c(0,0), c(-1,1))
Sigma <- matrix(c(1, 0.7, 0.7, 1), nr=2, nc=2) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus=mus, Sigmas=Sigmas, props=props, display="filled.contour", lwd=1)

## trivariate 
mus <- rbind(c(0,0,0), c(-1,0.5,1.5))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus=mus, Sigmas=Sigmas, props=props, dfs=c(11,8), dist="t")
}

\keyword{hplot}
