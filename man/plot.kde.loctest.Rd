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
  \item{x}{object of class \code{kde.loctest} (output from \code{\link{kde.local.test}})}
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{lcol}}{colour for KDE curve (1-d)}
      \item{\code{col}}{vector of 2 colours. First colour: sample 1>sample 2, second colour:
	sample 1<sample2.}
      \item{\code{add}}{flag to add to current plot. Default is FALSE.}
      \item{\code{rugsize}}{height of rug-like plot (1-d)}
      \item{\code{add.legend}}{flag to add legend. Default is FALSE.}
      \item{\code{pos.legend}}{position label for legend (1-d, 2-d)}
      \item{\code{add.contour}}{flag to add contour lines. Default is FALSE (2-d).} 
      \item{\code{alphavec}}{vector of transparency values for contour (3-d)}    
    }
    
    and those used in \code{\link{plot.kde}}
  }
} 

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is
  sent to graphics/RGL window. 
}

\details{
  For \code{kde.loctest} objects, the function headers are  
  \preformatted{
   ## univariate
   plot(x, lcol, col, add=FALSE, xlab="x", ylab, rugsize, add.legend=TRUE, 
        pos.legend="topright", alpha=1, ...)
   
   ## bivariate
   plot(x, col, add=FALSE, xlab="x", ylab="y", add.contour=FALSE, 
        add.legend=TRUE, pos.legend="topright", alpha=1, ...)

   ## trivariate 
   plot(x, col, color, add=FALSE, box=TRUE, axes=TRUE, alphavec=c(0.5, 0.5), 
        add.legend=TRUE, ...)}
}
 
\seealso{\code{\link{kde.local.test}}}

\examples{
## bivariate
data(air)
air.var <- c("co2","pm10","no")
air <- air[, c("date","time",air.var)]
air2 <- reshape(air, idvar="date", timevar="time", direction="wide")
a1 <- as.matrix(na.omit(air2[, paste0(air.var, ".08:00")]))
a2 <- as.matrix(na.omit(air2[, paste0(air.var, ".20:00")]))
colnames(a1) <- air.var
colnames(a2) <- air.var
loct <- kde.local.test(x1=a1[,c("co2","pm10")], x2=a2[,c("co2","pm10")])
plot(loct, lwd=1)

\donttest{## trivariate
loct <- kde.local.test(x1=a1, x2=a2)
plot(loct, xlim=c(0,800), ylim=c(0,300), zlim=c(0,300))}
}

\keyword{hplot}
