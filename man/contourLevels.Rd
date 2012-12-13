\name{contourLevels}
\alias{contourLevels}
\alias{contourLevels.kde}
\alias{contourLevels.kda.kde}
\alias{contourSizes}

\title{Contours for kde and kda.kde objects}
\description{
  Contour levels and sizes for \code{kde} and \code{kda.kde} objects.
} 

\usage{
contourLevels(x, ...)
\method{contourLevels}{kde}(x, prob, cont, nlevels=5, approx=FALSE, ...)
\method{contourLevels}{kda.kde}(x, prob, cont, nlevels=5, approx=FALSE, ...)

contourSizes(x, abs.cont, cont=c(25,50,75), approx=FALSE)
}

\arguments{
  \item{x}{an object of class \code{kde} or \code{kda.kde}}
  \item{prob}{vector of probabilities corresponding to highest density regions}
  \item{cont}{vector of percentages which correspond to the complement
	of \code{prob}}
  \item{abs.cont}{vector of absolute contour levels}
  \item{nlevels}{number of pretty contour levels} 
  \item{approx}{flag to compute approximate contour levels. Default is FALSE.} 
  \item{...}{other parameters}
}
  

\value{
  --For \code{contourLevels}, for \code{kde} objects, returns vector of heights.  For \code{kda.kde}
  objects, returns a list of vectors, one for each training group.

  --For \code{contourSizes}, an approximation of the Lebesgue measure of 
  level set, i.e. length (d=1), area (d=2), volume (d=3), hyper-volume (d>4). 
}



\details{
  The most straightforward is to specify \code{prob}. Heights of
  the corresponding highest density region with probability \code{prob} are
  computed. The \command{cont} parameter here is consistent with
  \command{cont} parameter from \command{plot.kde} and \command{plot.kda.kde}
  i.e. \code{cont=(1-prob)*100}\%. 
  If both \code{prob} and \code{cont} are missing then a pretty set of
  \code{nlevels} contours are computed.

  If \code{approx=FALSE}, then the exact KDE at \eqn{x}{x} is computed. Otherwise the exact 
  KDE is replaced by the weighted sum of the KDE at the nearest 
  grid points. This can dramatically reduce computation time for large data sets.
  
  The approximate Lebesgue measures are approximated by Riemann sums. Thsese are  rough approximations and depend highly on the estimation grid, and so should
  be interpreted carefully.
}
 

\seealso{\code{\link{contour}}, \code{\link{contourLines}}}

\examples{
## contourLevels kde 
x <- rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1)
Hx <- Hpi(x)
fhatx <- kde(x=x, H=Hx)
lev1 <- contourLevels(fhatx, prob=c(0.25, 0.5, 0.75))
lev2 <- contourLevels(fhatx, cont=c(75, 50, 25))      ## lev1==lev2

## contourLevels kda.kde
library(MASS)
data(iris)
ir <- iris[,1]
ir.gr <- iris[,5]
kda.fhat <- kda.kde(ir, ir.gr, hs=sqrt(c(0.01, 0.04, 0.07)))
contourLevels(kda.fhat, prob=c(0.25, 0.5, 0.75))

## contourSizes kde
x <- rmvnorm.mixt(n=1000, mus=c(0,0), Sigmas=diag(2), props=1)
Hx <- Hpi(x, binned=TRUE)
fhat <- kde(x=x, H=Hx, binned=TRUE)
plot(fhat, cont=25, xlim=c(-1,1), ylim=c(-1,1), asp=1) 
contourSizes(fhat, cont=25, approx=TRUE) 
 ## compare to approx circle of radius=0.75, vol=pi*0.75^2=1.767146
}

\keyword{hplot}
