\name{plot.kfs}
\alias{plot.kfs}
\title{Plot for kernel feature significance}

\description{
  Plot for kernel significant regions for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kfs}(x, display="filled.contour", col=7, colors, abs.cont,
   alpha=1, alphavec=0.4, add=FALSE, ...)
}

\arguments{
  \item{x}{object of class \code{kfs} (output from
    \code{\link{kfs}})}
  \item{display}{type of display, "slice" for contour plot,
    "persp" for perspective plot, "image" for image plot, "filled.contour"
    for filled contour plot (2-d); "plot3D", "rgl" (3-d)}
  \item{col,colors}{colour for contour region}
  \item{abs.cont}{absolute contour height. Default is 0.5.}
  \item{alpha}{transparency value for contour (2-d)}
  \item{alphavec}{vector of transparency values for contour (3-d)}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{...}{other graphics parameters used in \code{\link{plot.kde}}}
}
  
\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  graphics/RGL window.
}

\seealso{\code{\link{plot.kde}}}

\examples{
\donttest{library(MASS)
data(geyser)
geyser.fs <- kfs(geyser, binned=TRUE)
plot(geyser.fs)
}}

\keyword{hplot}
