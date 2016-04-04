\name{contour}
\alias{contourLevels}
\alias{contourLevels.kde}
\alias{contourLevels.kda}
\alias{contourSizes}

\title{Contours functions}
\description{
  Contour levels and sizes.
} 

\usage{
contourLevels(x, ...)
\method{contourLevels}{kde}(x, prob, cont, nlevels=5, approx=TRUE, ...)
\method{contourLevels}{kda}(x, prob, cont, nlevels=5, approx=TRUE, ...)

contourSizes(x, abs.cont, cont=c(25,50,75), approx=TRUE)
}

\arguments{
  \item{x}{an object of class \code{kde} or \code{kda}}
  \item{prob}{vector of probabilities corresponding to highest density regions}
  \item{cont}{vector of percentages which correspond to the complement
	of \code{prob}}
  \item{abs.cont}{vector of absolute contour levels}
  \item{nlevels}{number of pretty contour levels} 
  \item{approx}{flag to compute approximate contour levels. Default is TRUE.} 
  \item{...}{other parameters}
}
  

\value{
  --For \code{contourLevels}, for \code{kde} objects, returns vector of heights.  For \code{kda}
  objects, returns a list of vectors, one for each training group.

  --For \code{contourSizes}, an approximation of the Lebesgue measure of 
  level set, i.e. length (d=1), area (d=2), volume (d=3), hyper-volume (d>4). 
}



\details{
  --For \code{contourLevels}, the most straightforward is to specify \code{prob}. Heights of
  the corresponding highest density region with probability \code{prob} are
  computed. The \command{cont} parameter here is consistent with
  \command{cont} parameter from \command{plot.kde} and \command{plot.kda}
  i.e. \code{cont=(1-prob)*100}\%. 
  If both \code{prob} and \code{cont} are missing then a pretty set of
  \code{nlevels} contours are computed.
 
  --For \code{contourSizes}, the approximate Lebesgue measures are approximated by Riemann sums. Thsese are  rough approximations and depend highly on the estimation grid, and so should
  be interpreted carefully.

  If \code{approx=FALSE}, then the exact KDE is computed. Otherwise
  it is interpolated from an existing KDE grid. This can dramatically
  reduce computation time for large data sets. 
}
 

\seealso{\code{\link{contour}}, \code{\link{contourLines}}}

\examples{
set.seed(8192)
x <- rmvnorm.mixt(n=1000, mus=c(0,0), Sigmas=diag(2), props=1)
fhat <- kde(x=x, binned=TRUE)
contourLevels(fhat, cont=c(75, 50, 25))
contourSizes(fhat, cont=25, approx=TRUE) 
   ## compare to approx circle of radius=0.75 with area=1.77
}

\keyword{hplot}
