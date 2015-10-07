\name{plot.kcde}
\alias{plot.kcde}

\title{Plot for kernel cumulative distribution estimate}
\description{
  Plot for kernel cumulative distribution estimate 1- to 3-dimensional data.
}

\usage{
\method{plot}{kcde}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{kcde} (output from \code{\link{kcde}})}
  \item{...}{other graphics parameters used in \code{\link{plot.kde}}}
}
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  RGL window (not yet implemented). 
}

\details{
  For \code{kcde} objects, the function headers for the different dimensional data are  
  \preformatted{
  ## univariate
  plot(Fhat, xlab, ylab="Distribution function", add=FALSE, drawpoints=FALSE, 
       col.pt="blue", jitter=FALSE, ...) 

  ## bivariate
  plot(Fhat, display="persp", cont=seq(10,90, by=10), abs.cont, xlab, ylab, 
       zlab="Distribution function", cex=1, pch=1, add=FALSE, drawpoints=FALSE, 
       drawlabels=TRUE, theta=-30, phi=40, d=4, col.pt="blue", col, col.fun, 
       lwd=1, border=NA, thin=1, ...) 
  }
}
  
\seealso{\code{\link{plot.kde}}}

\examples{
library(MASS)
data(iris)
Fhat <- kcde(x=iris[,1])
plot(Fhat, xlab="Sepal.Length")
Fhat <- kcde(x=iris[,1:2])
plot(Fhat, thin=3)
}


\keyword{hplot}
