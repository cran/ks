\name{ksupp}
\alias{ksupp}

\title{Kernel support estimate}

\description{
  Kernel support estimate for 2-dimensional data.
}
\usage{
ksupp(fhat, cont=95, abs.cont, convex.hull=FALSE)
}

\arguments{
  \item{fhat}{object of class \code{kde}}
  \item{cont}{percentage for contour level curve. Default
  is 95.}
  \item{abs.cont}{absolute density estimate height for contour level curve}
  \item{convex.hull}{flag to compute convex hull of contour level
  curve. Default is FALSE.}
}
  

\value{
  A kernel support estimate is a 2-column matrix which delimits the
  (convex hull of the) level set of the density estimate \code{fhat}.  
}

\details{
  The kernel support estimate is the level set of the density estimate
  that exceeds the \code{cont} percent contour level. If this level set
  is a simply connected region, then this can suffice to be a
  conservative estimate of the density support. Otherwise, the convex
  hull of the level set is advised.  
}

\seealso{\code{\link{kde}}}

\examples{
library(oz)
data(grevillea)
fhat <- kde(x=grevillea)
fhat.supp <- ksupp(fhat, convex.hull=TRUE)
plot(fhat, display="filled.contour", cont=seq(10,90,by=10), drawlabels=FALSE)
plot(fhat, cont=95, add=TRUE)
polygon(fhat.supp, lty=2)
}

\keyword{smooth}
