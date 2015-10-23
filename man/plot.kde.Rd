\name{plot.kde}
\alias{plot.kde}

\title{Plot for kernel density estimate}
\description{
  Plot for kernel density estimate for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kde}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{kde} (output from \code{\link{kde}})}
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{display}}{type of display, "slice" for contour plot,
	"persp" for perspective plot, "image" for image plot, "filled.contour"
	for filled contour plot (1st form), "filled.contour2" (2nd form) (2-d)}
      \item{\code{cont}}{vector of percentages for contour level curves}
      \item{\code{abs.cont}}{vector of absolute density estimate heights for contour level curves}
      \item{\code{approx.cont}}{flag to compute approximate contour levels. Default is FALSE.}
      \item{\code{col}}{plotting colour for density estimate (1-d, 2-d)}
      \item{\code{col.cont}}{plotting colour for contours}
      \item{\code{col.fun}}{plotting colour function for contours}
      \item{\code{col.pt}}{plotting colour for data points}
      \item{\code{colors}}{vector of colours for each contour (3-d)}
      \item{\code{jitter}}{flag to jitter rug plot (1-d). Default is TRUE.}
      \item{\code{lwd.fc}}{line width for filled contours (2-d)}
      \item{\code{xlim,ylim,zlim}}{axes limits}
      \item{\code{xlab,ylab,zlab}}{axes labels}
      \item{\code{add}}{flag to add to current plot. Default is FALSE.}
      \item{\code{theta,phi,d,border}}{graphics parameters for perspective plots (2-d)}
      \item{\code{drawpoints}}{flag to draw data points on density estimate. Default is FALSE.}
      \item{\code{drawlabels}}{flag to draw contour labels (2-d). Default is TRUE.} 
      \item{\code{alpha}}{transparency value of plotting symbol (3-d)}
      \item{\code{alphavec}}{vector of transparency values for contours (3-d)}
      \item{\code{size}}{size of plotting symbol (3-d).}
    }
  }
}
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  RGL window. 
}


\details{
  For \code{kde} objects, the function headers for the different dimensional data are  
  \preformatted{
  ## univariate
  plot(fhat, xlab, ylab="Density function", add=FALSE, drawpoints=FALSE,
     col.pt="blue", col.cont=1, cont.lwd=1, jitter=FALSE, cont, abs.cont, 
     approx.cont=TRUE, ...) 
  
  ## bivariate
  plot(fhat, display="slice", cont=c(25,50,75), abs.cont, approx.cont=TRUE, 
  xlab, ylab, zlab="Density function", cex=1, pch=1, add=FALSE,
  drawpoints=FALSE, drawlabels=TRUE, theta=-30, phi=40, d=4, col.pt="blue", 
  col, col.fun, lwd=1, border=1, thin=3, lwd.fc=5, ...)

  ## trivariate
  plot(fhat, cont=c(25,50,75), abs.cont, approx.cont=FALSE, colors,
    add=FALSE, drawpoints=FALSE, alpha, alphavec, xlab, ylab, zlab, 
    size=3, col.pt="blue", ...)   
  }

  The 1-d plot is a standard plot of a 1-d curve. If
  \code{drawpoints=TRUE} then a rug plot is added. If \code{cont} is specified,
  the horizontal line on the x-axis indicates the \code{cont}\% highest 
  density level set.  
  
  There are different types of plotting displays for 2-d data available,
  controlled by the \code{display} parameter.
  (a) If \code{display="slice"} then a slice/contour plot
  is generated using \code{contour}.  
  (b) If \code{display} is \code{"filled.contour"} or \code{"filled.contour2"} 
  then a filled contour plot is generated.
  The default contours are at 25\%, 50\%, 75\% or
  \code{cont=c(25,50,75)} which are upper percentages of
  highest density regions.
  (c) If \code{display="persp"} then a perspective/wire-frame plot
  is generated.  The default z-axis limits \code{zlim} are the default
  from the usual \code{persp} command.
  (d) If \code{display="image"} then an image plot
  is generated. Default colours are the default from the usual
  \code{image} command.

  For 3-dimensional data, the interactive plot is a series of nested
  3-d contours. 
  The default contours are \code{cont=c(25,50,75)}.  The
  default \code{colors} are \code{heat.colors} and the
  default opacity \code{alphavec} ranges from 0.1 to 0.5.

  To specify contours, either one of \code{cont} or \code{abs.cont}
  is required. \code{cont} specifies upper percentages which
  correspond to probability contour regions. If \code{abs.cont} is set
  to particular values, then contours at these levels are drawn.
  This second option is useful for plotting
  multiple density estimates with common contour levels. See
  \code{\link{contourLevels}} for details on computing contour levels.   
  If \code{approx=FALSE}, then the exact KDE is computed. Otherwise
  it is interpolated from an existing KDE grid. This can dramatically
  reduce computation time for large data sets. 
}
 
\examples{
library(MASS)
data(iris)

## univariate example
fhat <- kde(x=iris[,2])
plot(fhat, cont=50, col.cont="blue", cont.lwd=2, xlab="Sepal length")

## bivariate example
fhat <- kde(x=iris[,2:3])
plot(fhat, display="filled.contour2", cont=seq(10,90,by=10))
plot(fhat, display="persp", thin=3, border=1, col="white")
\donttest{
## trivariate example
fhat <- kde(x=iris[,2:4])
plot(fhat, drawpoints=TRUE)
}}


\keyword{hplot}
