\name{plot.kda.kde}
\alias{plot.kda.kde}

\title{Kernel discriminant analysis plot
for 1- to 3-dimensional data}
\description{
  Kernel discriminant analysis plot for 1- to 3-dimensional data.}
}

\synopsis{
\method{plot}{kda.kde}(x, y, y.group, drawpoints=FALSE, ...)
}

\usage{
## univariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, xlim, ylim,
    xlab="x", ylab="Weighted density function", drawpoints=FALSE,
    col, ptcol, jitter=TRUE, ...)

## bivariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, cont=c(25,50,75),
    abs.cont, xlim, ylim, xlab, ylab, drawpoints=FALSE,
    drawlabels=TRUE, col, partcol, ptcol, ...)

## trivariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, cont=c(25,50,75),
   abs.cont, colors, alphavec, xlab, ylab, zlab, drawpoints=FALSE,
   size=3, ptcol="blue", ...)

}

\arguments{
  \item{x}{an object of class \code{kda.kde} (output from
    \code{\link{kda.kde}})}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{prior.prob}{vector of prior probabilities}  
  \item{cont}{vector of percentages for contour
    level curves}
  \item{abs.cont}{vector of absolute density estimate heights for contour
    level curves}
  \item{xlim,ylim}{axes limits}
  \item{xlab,ylab,zlab}{axes labels}
  \item{drawpoints}{if TRUE then draw data points}
  \item{drawlabels}{if TRUE then draw contour labels (2-d plot)}
  \item{jitter}{if TRUE then jitter rug plot (1-d plot)}
  \item{ptcol}{vector of colours for data points of each group}
  \item{partcol}{vector of colours for partition classes (1-d, 2-d plot)}
  \item{col}{vector of colours for density estimates (1-d, 2-d plot)}
  \item{colors}{vector of colours for contours of density estimates (3-d plot)}
  \item{alphavec}{vector of transparency values - one for each contour
    (3-d plot)}
  \item{size}{size of plotting symbol (3-d plot)}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of 1-d and 2-d density estimates for discriminant analysis is
  sent to graphics window. Plot for 3-d is sent to RGL window.
}

\details{
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
 
\references{  

  Bowman, A.W. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Clarendon Press. Oxford.
  
  Simonoff, J. S., (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.}
}


\seealso{\code{\link{kda.kde}}, \code{\link{kda}}}

\examples{
library(MASS)
data(iris)

## univariate example
ir <- iris[,1]
ir.gr <- iris[,5]
kda.fhat <- kda.kde(ir, ir.gr, hs=sqrt(c(0.01, 0.04, 0.07)))
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
H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H)
plot(kda.fhat)    
}

\keyword{ hplot}

