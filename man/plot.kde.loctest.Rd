\name{plot.kde.loctest}
\alias{plot.kde.loctest}

\title{Plot for kernel local significant difference regions}
\description{
  Plot for kernel local significant difference regions for 1- to 3-dimensional data. 
}

\usage{
\method{plot}{kde.loctest}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{kde.loctest} (output from \code{\link{kde.local.test}})}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d kernel density estimates are sent to graphics window. Plot
  for 3-d is sent to RGL window.
}


\details{
  For \code{kde.loctest} objects, the function headers are  
  \preformatted{
   ## univariate
   plot(x, col, add=FALSE, xlab="x", ylab, rugsize, add.legend=TRUE, 
     pos.legend="topright, ...)
   
   ## bivariate
   plot(x, lcol, col, add=FALSE, xlab="x", ylab="y", add.contour=FALSE, 
     add.legend=TRUE, pos.legend="topright", ...)

   ## trivariate 
   plot(x, col, add=FALSE, xlab="x", ylab="y", zlab="z", box=TRUE, 
     axes=TRUE, alphavec=c(0.5, 0.5), ...)
 }
  The arguments are
  \describe{
  \item{\code{x}}{object of class \code{kde.loctest}, output from \code{kde.local.test}}
  \item{\code{lcol}}{colour for KDE curve (1-d)}
  \item{\code{col}}{vector of 2 colours. First colour: sample 1 > sample 2, second colour: sample 1 < sample 2. Default is c("purple", "darkgreen").}
  \item{\code{add}}{flag to add to existing plot. Default is FALSE.}
  \item{\code{rugsize}}{height of rug-like plot (1-d)} 
  \item{\code{add.legend}}{flag to add legend. Default is FALSE (1-d, 2-d).}
  \item{\code{pos.legend}}{position label for legend (1-d, 2-d).}
  \item{\code{add.contour}}{flag to add contour lines. Default is FALSE (2-d).}
  \item{\code{box,axes}}{flags to add box and axes to plot (3-d).}
  \item{\code{alphavec}}{vector of 1 or 2 transparency values of difference regions (3-d).} 
  } 
 
  See \code{\link{plot.kde}} for more details. 
}
 

\seealso{\code{\link{kde}}, \code{\link{kde.local.test}}}

\examples{
## See examples in ? kde.local.test
}


\keyword{hplot}
