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
  \item{...}{other parameters} 
}
  

\value{
  A truncated kernel density (derivative) estimate inherits the same
  object class as the input estimate. 
}

\details{
  A simple truncation is performed on the kernel estimator. All the
  points in the estimation grid which are outside of the regions
  delimited by \code{boundary} are set to 0, and their probability
  mass is distributed proportionally to the density (derivative) value.
  %The \code{boundary} is typically obtained as the output from
  %\code{polypath} from the \pkg{mgcv} package which determines
  %automatically the interior/exterior of a set of polygon paths.
}

\seealso{\code{\link{kde}}, \code{\link{kdde}}}

\examples{
library(oz)
data(grevillea)
wa.coast <- ozRegion(section=1)
wa.polygon <- cbind(wa.coast$lines[[1]]$x, wa.coast$lines[[1]]$y)
fhat <- kde(x=grevillea)
fhat <- kde.truncate(fhat, wa.polygon)
oz(section=1, xlim=c(114.75, 121.5),  ylim=c(-33, -31.5))
plot(fhat, add=TRUE, cont=seq(10,90,by=10), col=2, drawlabels=FALSE, drawpoints=TRUE)
}

\keyword{smooth}
