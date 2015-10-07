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
  \item{...}{other graphics parameters:
    \describe{
      \item{\code{lcol}}{colour for KDE curve (1-d)}
      \item{\code{col}}{vector of 2 colours. Default is c("purple",
	"darkgreen"). First colour: sample 1>sample 2, second colour:
	sample 1<sample2.}
      \item{\code{add}}{flag to add to current plot. Default is FALSE.}
      \item{\code{rugsize}}{height of rug-like plot (1-d)}
      \item{\code{add.legend}}{flag to add legend. Default is FALSE (1-d, 2-d).}
      \item{\code{pos.legend}}{position label for legend (1-d, 2-d)}
      \item{\code{add.contour}}{flag to add contour lines. Default is FALSE (2-d).}     
    }
    
    and those used in \code{\link{plot.kde}}
  }
}
  
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is
  sent to RGL window. 
}



\details{
  For \code{kde.loctest} objects, the function headers are  
  \preformatted{
   ## univariate
   plot(x, lcol, col, add=FALSE, xlab="x", ylab, rugsize, add.legend=TRUE, 
     pos.legend="topright", ...)
   
   ## bivariate
   plot(x, col, add=FALSE, xlab="x", ylab="y", add.contour=FALSE, 
     add.legend=TRUE, pos.legend="topright", ...)

   ## trivariate 
   plot(x, col, add=FALSE, xlab="x", ylab="y", zlab="z", box=TRUE, axes=TRUE,
     alphavec=c(0.5, 0.5), ...)
   }
}
 

\seealso{\code{\link{kde.local.test}}}

\examples{
library(MASS)
data(crabs)
x1 <- crabs[crabs$sp=="B", c(4,6)]
x2 <- crabs[crabs$sp=="O", c(4,6)]
loct <- kde.local.test(x1=x1, x2=x2)
plot(loct)
}


\keyword{hplot}
