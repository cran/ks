\name{plotmixt}
\alias{plotmixt}

\title{Plot for 2- and 3-dimensional normal and t-mixture density functions.}
\description{
  Plot for 2- and 3-dimensional normal and t-mixture density functions.
}


\usage{
plotmixt(mus, Sigmas, props, dfs, dist="normal", ...)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{dfs}{vector of degrees of freedom}
  \item{dist}{"normal" - normal mixture, "t" - t-mixture}
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
plotmixt(mus, Sigmas, props)
plotmixt(mus, Sigmas, props, dfs=c(3,8), dist="t")

%%\dontrun{
## trivariate example
mus <- rbind(c(0,0,0), c(-1,1,1))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
plotmixt(mus, Sigmas, props)
plotmixt(mus, Sigmas, props,  dfs=c(3,8), dist="t")
%}
}

\keyword{ hplot}
