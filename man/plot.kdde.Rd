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
  \item{x}{an object of class \code{kdde} (output from \code{\link{kdde}})}
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{which.deriv.ind}}{index of the partial derivative to
	be plotted (>1-d)}
    }

    and those used in \code{\link{plot.kde}}}
}
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  RGL window.

  In addition to the display options inherited from \code{plot.kde}, the
  first derivative has \code{display="quiver"}. This is a quiver plot
  where the size and direction of the arrow indicates the
  magnitude/direction of the density gradient. See \code{quiver2D} from
  the \pkg{OceanView} package for more details. 
}


\details{
  For \code{kdde} objects, the function headers for the different dimensional data are  
  \preformatted{
  ## univariate
  plot(fhat, ylab="Density derivative function", ...)

  ## bivariate
  plot(fhat, which.deriv.ind=1, cont=c(25,50,75), abs.cont, display="slice", 
     zlab="Density derivative function", ...)
  }
}
 

\seealso{\code{\link{plot.kde}}}

\examples{
library(OceanView)
## univariate example
data(tempb)
fhat1 <- kdde(x=tempb[,"tmin"], deriv.order=1)   ## gradient [df/dx, df/dy]
plot(fhat1, xlab="Min. temp.")                   ## df/dx
points(20,predict(fhat1, x=20))

## bivariate example
fhat1 <- kdde(x=tempb[,c("tmin", "tmax")], deriv.order=1)
plot(fhat1, display="quiver")
  ## gradient [df/dx, df/dy]

\donttest{
fhat2 <- kdde(x=tempb[,c("tmin", "tmax")], deriv.order=2)
plot(fhat2, which.deriv.ind=2, display="persp", phi=15)
plot(fhat2, which.deriv.ind=2, display="filled.contour", col.fun=topo.colors)
  ## d^2 f/(dx dy): purple=-ve, green=zero, beige=+ve
s2 <- kcurv(fhat2)
plot(s2, display="filled.contour")
  ## summary curvature 
}}


\keyword{hplot}
