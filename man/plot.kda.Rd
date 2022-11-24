\name{plot.kda}
\alias{plot.kda}
\title{Plot for kernel discriminant analysis}

\description{
  Plot for kernel discriminant analysis for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kda}(x, y, y.group, ...)
}

\arguments{
  \item{x}{ object of class \code{kda} (output from \code{\link{kda}})}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{rugsize}}{height of rug-like plot for partition classes (1-d)}
      \item{\code{prior.prob}}{vector of prior probabilities}
      \item{\code{col.part}}{vector of colours for partition classes (1-d, 2-d)}
    }
    
    and those used in \code{\link{plot.kde}}
  }
}
  
\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  graphics/RGL window. 
}

\details{
  For \code{kda} objects, the function headers for the different dimensional data are 
  \preformatted{
  ## univariate
  plot(x, y, y.group, prior.prob=NULL, xlim, ylim, xlab, 
       ylab="Weighted density function", drawpoints=FALSE, col, col.fun, 
       col.part, col.pt, lty, jitter=TRUE, rugsize, add=FALSE, alpha=1, ...)

  ## bivariate
  plot(x, y, y.group, prior.prob=NULL, display.part="filled.contour",
       cont=c(25,50,75), abs.cont, approx.cont=TRUE, xlim, ylim, xlab, ylab,
       drawpoints=FALSE, drawlabels=TRUE, cex=1, pch, lty, part=TRUE, col, 
       col.fun, col.part, col.pt, alpha=1, lwd=1, lwd.part=0, add=FALSE, ...)

  ## trivariate
  plot(x, y, y.group, prior.prob=NULL, display="plot3D", cont=c(25,50,75), 
       abs.cont, approx.cont=TRUE, colors, col, col.fun, col.pt, alpha=0.5, 
       alphavec, xlab, ylab, zlab, drawpoints=FALSE, size=3, cex=1, pch, 
       theta=-30, phi=40, d=4, ticktype="detailed", bty="f", add=FALSE, ...)}
}
 
\seealso{\code{\link{kda}}, \code{\link{kde}}}

\examples{
library(MASS)
data(iris)

## univariate example
ir <- iris[,1]
ir.gr <- iris[,5]
kda.fhat <- kda(x=ir, x.group=ir.gr, xmin=3, xmax=9)
plot(kda.fhat, xlab="Sepal length")

## bivariate example
ir <- iris[,1:2]
ir.gr <- iris[,5]
kda.fhat <- kda(x=ir, x.group=ir.gr)
plot(kda.fhat, alpha=0.2, drawlabels=FALSE)

## trivariate example
ir <- iris[,1:3]
ir.gr <- iris[,5]
kda.fhat <- kda(x=ir, x.group=ir.gr)
plot(kda.fhat) 
   ## colour=species, transparency=density heights}

\keyword{hplot}
