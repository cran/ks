\name{plot.kda.kde}
\alias{plot.kda.kde}

\title{Kernel discriminant analysis plot
for 1- to 3-dimensional data}
\description{
  Kernel discriminant analysis plot for 1- to 3-dimensional data.}
}

\synopsis{
\method{plot}{kda.kde}(x, y, y.group, ...)
}

\usage{
## univariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, xlim, ylim,
    xlab="x", ylab="Weighted density function", drawpoints=TRUE,
    lty, lcol, col, ptcol, ...)

## bivariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, cont=c(25,50,75),
    abs.cont, xlim, ylim, xlab, ylab, drawpoints=FALSE,
    drawlabels=TRUE, cex=1, pch, lty, col, lcol, ptcol, ...)

## trivariate
\method{plot}{kda.kde}(x, y, y.group, prior.prob=NULL, cont=c(25,50), colors,
   alphavec, xlab, ylab, zlab, drawpoints=FALSE, size=3,
   ptcol="blue", ...)

}

\arguments{
  \item{x}{an object of class \code{kda.kde} (output from
    \code{\link{kda.kde}})}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{prior.prob}{vector of prior probabilities}  
  \item{cont}{vector of percentiles of density estimate heights for contour
    level curves}
  \item{abs.cont}{vector of absolute density estimate heights for contour
    level curves - only one of \code{cont} or \code{abs.cont} needs to
	be specified}
  \item{cex,pch,lty,xlim,ylim,xlab,ylab,zlab}{usual graphics parameters}
  \item{drawpoints}{if TRUE then draw data points}
  \item{drawlabels}{if TRUE then draw contour labels (2-d plot)}
  \item{col}{vector of colours for partition classes}
  \item{ptcol}{vector of colours for data points of each group}
  \item{lcol}{vector of colours for contour lines of density estimates}
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
  If \code{prior.prob} is set to a particular value then this is used. 
  The default is \code{NULL} which means that the sample proportions are used.

  The object \code{x} contains the training data and its group
  labels.  If \code{y} and \code{y.group} are missing then the training
  data points are plotted.  Otherwise, the test data \code{y} are plotted.

  -- For 1-d plots:
  
  The partition induced by the discriminant analysis is plotted as rug
  plot (with the ticks inside the axes). If \code{drawpoints=TRUE} then
  the data points are plotted as a rug plot with the ticks outside the
  axes, their colour is controlled by \code{ptcol}.
  
  -- For 2-d plots:
  
  If \code{display="part"} then a partition induced by the discriminant
  analysis is also plotted. If this is not desired, set
  \code{display=""}. Its colours are controlled by \code{col}
  (the default is \code{heat.colors}).
  The plotting symbols are set by \code{pch} and the colour by \code{ptcol}.
  The default contours are at 25\%, 50\%, 75\% or
  \code{cont=c(25,50,75)} for the highest
  density regions, see Hyndman (1996).
  The line types are set by \code{lty}. 
  

  -- For 3-d plots:
  
  Default contours are \code{cont=c(25,50)} for the highest density regions,
  see Hyndman (1996). The colour for each group
  is set using \code{colors} - default is \code{heat.colors}. The
  transparency of 
  each contour (within each group) is \code{alphavec}, with the default
  range from 0.1 to 0.5. 
}
 
\references{  

  Bowman, A.W. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Clarendon Press. Oxford.

  Hyndman, R.J. (1996) Computing and graphing highest density regions,
    \emph{American Statistician}, \bold{50}, 120--126

  Simonoff, J. S., (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.}
}


\seealso{\code{\link{kda.kde}}, \code{\link{kda}}}

\examples{
library(MASS)
data(iris)

xlim <- c(4,8)
ylim <- c(2,4.5)

## univariate example
ir <- iris[,1]
ir.gr <- iris[,5]

kda.fhat <- kda.kde(ir, ir.gr, hs=sqrt(c(0.01, 0.04, 0.07)))
plot(kda.fhat, xlab="Sepal length", ptcol=1:3)

## bivariate example
ir <- iris[,1:2]
ir.gr <- iris[,5]
H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H)

plot(kda.fhat, cont=0)
plot(kda.fhat, xlim=xlim, ylim=ylim, col=c("transparent", "grey", "#8f8f8f"))

## trivariate example

ir <- iris[,1:3]
ir.gr <- iris[,5] 
H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.fhat <- kda.kde(ir, ir.gr, Hs=H)
plot(kda.fhat)
}
\keyword{ hplot}

