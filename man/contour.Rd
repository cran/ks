\name{contour}
\alias{contourLevels}
\alias{contourLevels.kde}
\alias{contourLevels.kda}
\alias{contourLevels.kdde}
\alias{contourSizes}
\alias{contourProbs}
\title{Contour functions}

\description{
  Contour levels and sizes.
} 

\usage{
contourLevels(x, ...)
\method{contourLevels}{kde} (x, prob, cont, nlevels=5, approx=TRUE, ...)
\method{contourLevels}{kda} (x, prob, cont, nlevels=5, approx=TRUE, ...)
\method{contourLevels}{kdde}(x, prob, cont, nlevels=5, approx=TRUE, which.deriv.ind=1, ...) 

contourSizes(x, abs.cont, cont=c(25,50,75), approx=TRUE)
contourProbs(x, abs.cont, cont=c(25,50,75), approx=TRUE)
}

\arguments{
  \item{x}{object of class \code{kde}, \code{kdde} or \code{kda}}
  \item{prob}{vector of probabilities corresponding to highest density regions}
  \item{cont}{vector of percentages which correspond to the complement
	of \code{prob}}
  \item{abs.cont}{vector of absolute contour levels}
  \item{nlevels}{number of pretty contour levels} 
  \item{approx}{flag to compute approximate contour levels. Default is
    TRUE.} 
  \item{which.deriv.ind}{partial derivative index. Default is 1.} 
  \item{...}{other parameters}
}
  
\value{
  --For \code{contourLevels}, for \code{kde} objects, returns vector of
  heights.  For \code{kda} objects, returns a list of vectors, one for
  each training group. For \code{kdde} objects, returns a matrix of
  vectors, one row for each partial derivative. 

  --For \code{contourSizes}, returns an approximation of the Lebesgue measure of 
  level set, i.e. length (d=1), area (d=2), volume (d=3), hyper-volume (d>4).
  
  --For \code{contourProbs}, returns an approximation of the probability measure of 
  level set. 
}

\details{
  --For \code{contourLevels}, the most straightforward is to specify \code{prob}.
  The heights of the corresponding highest density region with probability \code{prob} are
  computed. The \command{cont} parameter here is consistent with
  \command{cont} parameter from \command{plot.kde}, \command{plot.kdde}, and \command{plot.kda}
  i.e. \code{cont=(1-prob)*100\%}. 
  If both \code{prob} and \code{cont} are missing then a pretty set of
  \code{nlevels} contours are computed.
 
  --For \code{contourSizes}, the length, area, volume etc. and for \code{contourProbs}, 
  the probability, are approximated by Riemann sums. These are rough approximations and
  depend highly on the estimation grid, and so should be interpreted carefully. 

  If \code{approx=FALSE}, then the exact KDE is computed. Otherwise
  it is interpolated from an existing KDE grid: this can dramatically
  reduce computation time for large data sets. 
}
 
\seealso{\code{\link{contour}}, \code{\link{contourLines}}}

\examples{
set.seed(8192)
x <- rmvnorm.mixt(n=1000, mus=c(0,0), Sigmas=diag(2), props=1)
fhat <- kde(x=x, binned=TRUE)
contourLevels(fhat, cont=c(75, 50, 25))
contourProbs(fhat, abs.cont=contourLevels(fhat, cont=50))
  ## compare approx prob with target prob=0.5
contourSizes(fhat, cont=25, approx=TRUE) 
   ## compare to approx circle of radius=0.75 with area=1.77
}

\keyword{hplot}
