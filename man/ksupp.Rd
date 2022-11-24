\name{ksupp}
\alias{ksupp}
\alias{plot.ksupp}
\title{Kernel support estimate}

\description{
  Kernel support estimate for 2 and 3-dimensional data.
}
\usage{
ksupp(fhat, cont=95, abs.cont, convex.hull=TRUE)

\method{plot}{ksupp}(x, display="plot3D", ...)
}

\arguments{
  \item{fhat}{object of class \code{kde}}
  \item{cont}{percentage for contour level curve. Default
  is 95.}
  \item{abs.cont}{absolute density estimate height for contour level curve}
  \item{convex.hull}{flag to compute convex hull of contour level
  curve. Default is TRUE.}
  \item{x}{object of class \code{ksupp}}
  \item{display}{one of "plot3D", "rgl" (required for 3-d only)}
  \item{...}{other graphics parameters}
}
  
\value{
  A kernel support estimate is an object of class \code{ksupp}, i.e. a 2- or 3-column matrix which delimits the (convex hull of the) level set of the density estimate \code{fhat}.  
}

\details{
  The kernel support estimate is the level set of the density estimate
  that exceeds the \code{cont} percent contour level. If this level set
  is a simply connected region, then this can suffice to be a
  conservative estimate of the density support. Otherwise, the convex
  hull of the level set is advised. For 2-d data, the convex hull is computed by \code{chull}; for 3-d data, it is computed by \code{geometry::convhulln}.
}

\seealso{\code{\link{kde}}}

\examples{
library(oz)
data(grevillea)
fhat <- kde(x=grevillea)
fhat.supp <- ksupp(fhat)
plot(fhat, display="filled.contour", cont=seq(10,90,by=10))
plot(fhat, cont=95, add=TRUE, col=1)
plot(fhat.supp,  lty=2)

data(iris)
fhat <- kde(x=iris[,1:3])
fhat.supp <- ksupp(fhat)
plot(fhat)
plot(fhat.supp, add=TRUE, col=3, alpha=0.1)
}

\keyword{smooth}
