\name{plot.kda.kde}
\alias{plot.kda.kde}

\title{Plot for kernel discriminant analysis}
\description{
  Plot for kernel discriminant analysis for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kda.kde}(x, y, y.group, ...)

}

\arguments{
  \item{x}{ object of class \code{kda.kde} (output from \code{\link{kda.kde}})}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d kernel density estimates are sent to graphics window. Plot
  for 3-d is sent to RGL window.
}


\details{
  For \code{kda.kde} objects, the function headers for the different dimensional data are 
  \preformatted{
  ## univariate
  plot(x, y, y.group, prior.prob=NULL, xlim, ylim, xlab="x", 
    ylab="Weighted density function", drawpoints=FALSE, col, ptcol, 
    jitter=TRUE, rugsize, ...)

  ## bivariate
  plot(x, y, y.group, prior.prob=NULL, cont=c(25,50,75), abs.cont, 
    approx.cont=FALSE, xlim, ylim, xlab, ylab, drawpoints=FALSE, 
    drawlabels=TRUE, col, partcol, ptcol, ...)

  ## trivariate
  plot(x, y, y.group, prior.prob=NULL, cont=c(25,50,75), abs.cont, 
    approx.cont=FALSE, colors, alphavec, xlab, ylab, zlab, 
    drawpoints=FALSE, size=3, ptcol="blue", ...)
  }

 
  The arguments are
  \describe{
  \item{\code{cont}}{vector of percentages for contour level curves}
  \item{\code{abs.cont} }{ vector of absolute density estimate heights for contour level curves}
  \item{\code{approx.cont} }{ flag to compute approximate contour levels. Default is FALSE.} 
  \item{\code{ptcol}}{ plotting colour for data points}
  \item{\code{col}}{ plotting colour for density estimate (1-d, 2-d plot)}
  \item{\code{colors}}{ vector of colours for each contour (3-d plot)}
  \item{\code{jitter}}{ flag to jitter rug plot (1-d plot). Default is TRUE.}
  \item{\code{rugsize}}{height of rug-like plot for partition classes (1-d plot)} 
  \item{\code{xlim, ylim}}{ axes limits}
  \item{\code{xlab,ylab,zlab}}{ axes labels}
  \item{\code{add}}{ flag to add to current plot. Default is FALSE.}
  \item{\code{drawpoints}}{ flag to draw data points on density estimate. Default is FALSE.}
  \item{\code{drawlabels}}{ flag to draw contour labels (2-d plot). Default is TRUE.}
  \item{\code{alphavec}}{ vector of transparency values for contours (3-d plot)}
  \item{\code{size}}{ size of plotting symbol (3-d plot)}
  \item{\code{prior.prob} }{ vector of prior probabilities}  
  \item{\code{partcol} }{ vector of colours for partition classes (1-d, 2-d plot)} 
  }
  See \code{\link{plot.kde}} for more details.
}
 

\seealso{\code{\link{kda.kde}}, \code{\link{kde}}}

\examples{
library(MASS)
data(iris)

## univariate example
ir <- iris[,1]
ir.gr <- iris[,5]
hs <- hkda(x=ir, x.gr=ir.gr)
kda.fhat <- kda.kde(ir, ir.gr, hs=hs, xmin=3, xmax=9)
plot(kda.fhat, xlab="Sepal length")

## bivariate example
ir <- iris[,1:2]
ir.gr <- iris[,5]
H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H)
plot(kda.fhat)   
}


\keyword{hplot}
