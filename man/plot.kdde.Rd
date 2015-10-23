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
## univariate example
data(unicef)
fhat1 <- kdde(x=unicef[,1], deriv.order=1)  ## gradient [df/dx, df/dy]
plot(fhat1, xlab="Under-5")                 ## df/dx
points(30,predict(fhat1, x=30))

## bivariate example
fhat2 <- kdde(x=unicef, deriv.order=2)
plot(fhat2, which.deriv.ind=2, display="persp", phi=20)
plot(fhat2, which.deriv.ind=2, display="filled.contour2", col.fun=topo.colors)
  ## d^2 f/(dx dy): purple=-ve, green=zero, beige=+ve
}


\keyword{hplot}
