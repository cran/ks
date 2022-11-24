\name{plot.kdde}
\alias{plot.kdde}
\title{Plot for kernel density derivative estimate}

\description{
  Plot for kernel density derivative estimate for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kdde}(x, ...)
}

\arguments{
  \item{x}{object of class \code{kdde} (output from \code{\link{kdde}})}
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{which.deriv.ind}}{index of the partial derivative to
	be plotted (>1-d)}
    }

    and those used in \code{\link{plot.kde}}}
}
  
\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  graphics/RGL window.

  In addition to the display options inherited from \code{plot.kde}, the
  first derivative has \code{display="quiver"}. This is a quiver plot
  where the size and direction of the arrow indicates the
  magnitude/direction of the density gradient. See \code{quiver} from
  the \pkg{pracma} package for more details. 
}

\details{
  For \code{kdde} objects, the function headers for the different dimensional data are  
  \preformatted{
  ## univariate
  plot(fhat, ylab="Density derivative function", cont=50, abs.cont, alpha=1, ...)

  ## bivariate
  plot(fhat, which.deriv.ind=1, cont=c(25,50,75), abs.cont, display="slice", 
       zlab="Density derivative function", col, col.fun, alpha=1, kdde.flag=TRUE, 
       thin=3, transf=1, neg.grad=FALSE, ...)
  
  ## trivariate 
  plot(fhat, which.deriv.ind=1, display="plot3D", cont=c(25,50,75), abs.cont, 
       colors, col, col.fun, ...)}
}
 
\seealso{\code{\link{plot.kde}}}

\examples{
## univariate example
data(tempb)
fhat1 <- kdde(x=tempb[,"tmin"], deriv.order=1)   ## gradient [df/dx, df/dy]
plot(fhat1, xlab="Min. temp.", col.cont=4)       ## df/dx
points(20,predict(fhat1, x=20))

## bivariate example
fhat1 <- kdde(x=tempb[,c("tmin", "tmax")], deriv.order=1)
plot(fhat1, display="quiver")
  ## gradient [df/dx, df/dy]

\donttest{fhat2 <- kdde(x=tempb[,c("tmin", "tmax")], deriv.order=2)
plot(fhat2, which.deriv.ind=2, display="persp", phi=10)
plot(fhat2, which.deriv.ind=2, display="filled.contour")
  ## d^2 f/(dx dy): blue=-ve, red=+ve
s2 <- kcurv(fhat2)
plot(s2, display="filled.contour", alpha=0.5, lwd=1)
  ## summary curvature 

## trivariate example  
data(iris)
fhat1 <- kdde(iris[,2:4], deriv.order=1)
plot(fhat1)}}

\keyword{hplot}
