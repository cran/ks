\name{contourLevels}
\alias{contourLevels}
\alias{contourLevels.kde}
\alias{contourLevels.kda.kde}

\title{Contour levels for kde and kda.kde objects}
\description{
  Contour levels for \code{kde} and \code{kda.kde} objects.
} 

\usage{
contourLevels(x, ...)

\method{contourLevels}{kde}(x, prob, cont, nlevels=5, approx=FALSE, ...)

\method{contourLevels}{kda.kde}(x, prob, cont, nlevels=5, approx=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{kde} or \code{kda.kde}}
  \item{prob}{vector of probabilities corresponding to highest density regions}
  \item{cont}{vector of percentages which correspond to the complement
	of \code{prob}}
  \item{nlevels}{number of pretty contour levels} 
  \item{approx}{flag to compute approximate contour levels. Default is FALSE.} 
  \item{...}{other parameters}
}
  

\value{
  For \code{kde} objects, returns vector of heights.  For \code{kda.kde}
  objects, returns a list of vectors, one for each training group.
}



\details{
  The most straightforward is to specify \code{prob}. Heights of
  the corresponding highest density region with probability \code{prob} are
  computed. The \command{cont} parameter here is consistent with
  \command{cont} parameter from \command{plot.kde} and \command{plot.kda.kde}
  i.e. \code{cont=(1-prob)*100}\%. If both \code{prob} and \code{cont} are missing then a pretty set of
  \code{nlevels} contours are computed.

  If \code{approx=FALSE}, then the exact KDE at \eqn{x}{x} is computed. Otherwise the exact 
  KDE is replaced by the weighted sum of the KDE at the nearest 
  grid points. This can dramatically reduce computation time for large data sets.
}
 

\seealso{\code{\link{contour}}, \code{\link{contourLines}}}

\examples{
## kde 
x <- rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1)
Hx <- Hpi(x)
fhatx <- kde(x=x, H=Hx)
lev1 <- contourLevels(fhatx, prob=c(0.25, 0.5, 0.75))
lev2 <- contourLevels(fhatx, cont=c(75, 50, 25))      ## lev1==lev2

## kda.kde
library(MASS)
data(iris)
ir <- iris[,1]
ir.gr <- iris[,5]
kda.fhat <- kda.kde(ir, ir.gr, hs=sqrt(c(0.01, 0.04, 0.07)))
contourLevels(kda.fhat, prob=c(0.25, 0.5, 0.75))
}

\keyword{hplot}
