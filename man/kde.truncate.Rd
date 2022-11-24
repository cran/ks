\name{kde.truncate}
\alias{kde.truncate}
\alias{kdde.truncate}
\title{Truncated kernel density derivative estimate}

\description{
  Truncated kernel density derivative estimate for 2-dimensional data.
}

\usage{
kde.truncate(fhat, boundary) 
kdde.truncate(fhat, boundary) 
}

\arguments{
  \item{fhat}{object of class \code{kde} or \code{kdde}}
  \item{boundary}{two column matrix delimiting the boundary for truncation}
}
  
\value{
  A truncated kernel density (derivative) estimate inherits the same
  object class as the input estimate. 
}

\details{
  A simple truncation is performed on the kernel estimator. All the
  points in the estimation grid which are outside of the regions
  delimited by \code{boundary} are set to 0, and their probability
  mass is distributed proportionally to the remaining density (derivative) values.
}

\seealso{\code{\link{kde}}, \code{\link{kdde}}}

\examples{
data(worldbank)
wb <- as.matrix(na.omit(worldbank[,c("internet", "ag.value")]))
fhat <- kde(x=wb)
rectb <- cbind(x=c(0,100,100,0,0), y=c(0,0,100,100,0))
fhat.b <- kde.truncate(fhat, boundary=rectb)
plot(fhat, col=1, xlim=c(0,100), ylim=c(0,100))
plot(fhat.b, add=TRUE, col=4)
rect(0,0,100,100, lty=2)

library(oz)
data(grevillea)
wa.coast <- ozRegion(section=1)
wa.polygon <- cbind(wa.coast$lines[[1]]$x, wa.coast$lines[[1]]$y)
fhat1 <- kdde(x=grevillea, deriv.order=1)
fhat1 <- kdde.truncate(fhat1, wa.polygon)
oz(section=1, xlim=c(113,122), ylim=c(-36,-29))
plot(fhat1, add=TRUE, display="filled.contour")
}

\keyword{smooth}
