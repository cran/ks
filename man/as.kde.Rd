\name{as.kde}
\alias{as.kde}
\title{Quasi kernel density estimate}

\description{
  Quasi kernel density estimate for gridded data. 
}

\usage{
as.kde(x, density)
}

\arguments{
  \item{x}{tidy matrix with first d columns are evaluation points of an estimation grid, and last column is estimate value; or list with elements \code{eval.points}, \code{estimate}}
  \item{density}{flag to compute bona fide density. Default is TRUE.} 
}

\value{
  A quasi-kernel density estimate is an object of class \code{kqde} (for \code{density=TRUE}) or \code{kqdde} (for \code{density=FALSE}) which is a
  list with fields: 
  \item{x}{zero vector or matrix}
  \item{eval.points}{vector or list of points at which the estimate is evaluated - copied from input \code{x}} 
  \item{estimate}{density estimate - copied from input \code{x}}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{cont}{vector of probability contour levels}
  This mimics the output from \code{\link{kde}} (for \code{density=TRUE}) or \code{\link{kdde}}(for \code{density=FALSE}).   
}

\details{
  This function converts the input \code{x} into a quasi-kernel density estimate which is treated similarly to the output from \code{\link{kde}} or \code{\link{kdde}}. If \code{density=TRUE} then all the estimate values are non-negative. If \code{density=FALSE} then the estimate values can be positive or negative. A \code{contourLevels} method computes an approximation of the density contour levels, and a \code{plot} method computes analogous plots. 
}

\seealso{\code{\link{kde}}, \code{\link{kdde}}}

\examples{
## univariate 
ep1 <- seq(0,1, length=401)
est1 <- dbeta(ep1, shape1=2, shape2=5)
f1 <- data.frame(ep1, estimate=est1)
f1.kde <- as.kde(f1)
plot(f1.kde)

## bivariate
ep2 <- expand.grid(x=seq(-pi,pi, length=151), y=seq(-pi,pi, length=151))
est2 <- cos(ep2[,1]*pi/2) + sin(ep2[,2]*pi/2)  
f2 <- data.frame(ep2, estimate=est2)
f2.kde <- as.kde(f2, density=FALSE)
plot(f2.kde, display="filled.contour")
plot(f2.kde, display="persp", phi=10)

\donttest{## trivariate
mus3 <- c(0,0,0)
Sigmas3 <- 1/2*diag(c(1,1/3,2))
props3 <- 1
ep3 <- expand.grid(x=seq(-3,3, length=51), y=seq(-3,3, length=51), z=seq(-3,3, length=51))
est3 <- dmvnorm.mixt(ep3, mus=mus3, Sigmas=Sigmas3, props=props3)
f3 <- data.frame(ep3, estimate=est3)
f3.kde <- as.kde(f3)
plot(f3.kde)}
}

\keyword{smooth}
