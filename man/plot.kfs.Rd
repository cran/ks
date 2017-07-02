\name{plot.kfs}
\alias{plot.kfs}

\title{Plot for kernel feature significance}
\description{
  Plot for kernel significant regions for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kfs}(x, display="filled.contour2", col="orange", colors="orange", abs.cont,
   alphavec=0.4, add=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{kfs} (output from
    \code{\link{kfs}})}
  \item{display}{type of display, "slice" for contour plot,
    "persp" for perspective plot, "image" for image plot, "filled.contour"
    for filled contour plot (1st form), "filled.contour2" (2nd form) (2-d)}
  \item{col,colors}{colour for contour region (1-d, 2-d), (3-d)}
  \item{abs.cont}{absolute contour height. Default is 0.5.}
  \item{alphavec}{transparency value for contour (3-d)}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{...}{other graphics parameters used in \code{\link{plot.kde}}}
}
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  RGL window.
}

\seealso{\code{\link{plot.kde}}}

\examples{
library(MASS)
data(geyser)
geyser.fs <- kfs(geyser, binned=TRUE)
plot(geyser.fs)
}


\keyword{hplot}
