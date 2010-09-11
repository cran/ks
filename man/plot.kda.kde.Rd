\name{plot.kda.kde}
\alias{plot}
\alias{plot.kda.kde}

\title{Kernel discriminant analysis plot for 1- to 3-dimensional data}
\description{
  Kernel discriminant analysis plot for 1- to 3-dimensional data.
}

\usage{
\method{plot}{kda.kde}(x, y, y.group, ...)
}

\arguments{
  \item{x}{an object of class \code{kda.kde} (output from
    \code{\link{kda.kde}})}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d density estimates for discriminant analysis is
  sent to graphics window. Plot for 3-d is sent to RGL window.
}

\details{
  Function headers for the different dimensional data are
  \preformatted{
  ## univariate
  plot(x, y, y.group, prior.prob=NULL, xlim, ylim,
    xlab="x", ylab="Weighted density function", drawpoints=FALSE,
    col, ptcol, jitter=TRUE, ...)

  ## bivariate
  plot(x, y, y.group, prior.prob=NULL, cont=c(25,50,75),
    abs.cont, approx.cont=FALSE, xlim, ylim, xlab, ylab, drawpoints=FALSE,
    drawlabels=TRUE, col, partcol, ptcol, ...)

  ## trivariate
  plot(x, y, y.group, prior.prob=NULL, cont=c(25,50,75),
    abs.cont,  approx.cont=FALSE, colors, alphavec, xlab, ylab, zlab, 
    drawpoints=FALSE, size=3, ptcol="blue", ...)
  }
  The arguments are 
  \tabular{ll}{
  \code{prior.prob} \tab vector of prior probabilities \cr  
  \code{cont} \tab vector of percentages for contour level curves \cr
  \code{abs.cont} \tab vector of absolute density estimate heights for contour level curves \cr
  \code{approx.cont} \tab flag to compute approximate contour levels \cr
  \code{xlim, ylim} \tab axes limits \cr
  \code{xlab, ylab, zlab} \tab axes labels \cr
  \code{drawpoints} \tab flag to draw data points. Default is FALSE. \cr
  \code{drawlabels} \tab flag to draw contour labels (2-d plot). Default is TRUE.\cr
  \code{jitter} \tab flag to jitter rug plot (1-d plot). Default is TRUE. \cr
  \code{ptcol} \tab vector of colours for data points of each group \cr
  \code{partcol} \tab vector of colours for partition classes (1-d, 2-d plot) \cr
  \code{col} \tab vector of colours for density estimates (1-d, 2-d plot) \cr
  \code{colors} \tab vector of colours for contours of density estimates (3-d plot) \cr
  \code{alphavec} \tab vector of transparency values - one for each contour (3-d plot) \cr
  \code{size} \tab size of plotting symbol (3-d plot)
  }
  -- For 1-d plots:
  
  The partition induced by the discriminant analysis is plotted as rug
  plot (with the ticks inside the axes). If \code{drawpoints=TRUE} then
  the data points are plotted as a rug plot with the ticks outside the
  axes, their colour is controlled by \code{ptcol}.
  
  -- For 2-d plots:
  
  The partition classes are displayed using the colours in \code{partcol}.
  The default contours of the density estimate are  25\%, 50\%, 75\% or
  \code{cont=c(25,50,75)} for highest density regions.
  See \command{\link{plot.kde}} for more details.
  
  -- For 3-d plots:
  
  Default contours are \code{cont=c(25,50,75)} for highest density
  regions.  See \command{\link{plot.kde}} for more 
  details. The colour of each group is \code{colors}. The transparency of
  each contour (within each group) is \code{alphavec}. Default range is
  0.1 to 0.5.

  -- If \code{prior.prob} is set to a particular value then this is used. 
  The default is \code{NULL} which means that the sample proportions are used.

  If \code{y} and \code{y.group} are missing then the training
  data points are plotted.  Otherwise, the test data \code{y} are plotted.

}
 

\seealso{\code{\link{kda.kde}}, \code{\link{kda}}}

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
plot(kda.fhat, cont=0, partcol=4:6)
plot(kda.fhat, drawlabels=FALSE, drawpoints=TRUE)

## trivariate example
## colour indicates species, transparency indicates density heights
ir <- iris[,1:3]
ir.gr <- iris[,5] 
H <- Hkda(ir, ir.gr, bw="plugin")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H, compute.cont=TRUE)
plot(kda.fhat)   
}

\keyword{ hplot}

