\name{plot.dade}
\alias{plot.dade}

\title{Density estimates and partition plot for discriminant analysis
for bi-and tri-variate data}
\description{
  Density estimate plot and partition for bi- and tri-variate data for kernel, linear and
  quadratic discriminant analysis}
}

\synopsis{
plot.dade(x, y, y.group, prior.prob=NULL, display="part",
          cont=c(25,50,75), ncont=NULL, ...)
}

\usage{
## bivariate
\method{plot}{dade}(x, y, y.group, prior.prob=NULL, display="part",
    cont=c(25,50,75), ncont=NULL, xlim, ylim, xlabs="x", ylabs="y",
    drawlabels=TRUE, cex=1, pch, lty, col, lcol, ...)

## trivariate
\method{plot}{dade}(x, y, y.group, prior.prob=NULL, 
    cont=c(25,50), colors, alphalo=0.2, alphahi=0.6, ...)
}

\arguments{
  \item{x}{an object of class \code{dade} i.e. output from
    \code{\link{kda.kde}} or \code{\link{pda.pde}}}
  \item{display}{include plot of partition classes}
  \item{y}{matrix of test data points}
  \item{y.group}{vector of group labels for test data points}
  \item{prior.prob}{vector of prior probabilities}  
  \item{cont}{vector of percentages (of maximum height) for contour
    level curves}
  \item{ncont}{number of contour level curves}
  \item{cex,pch,lty,xlim, ylim, xlabs, ylabs}{usual graphics parameters}
  \item{drawlabels}{draw contour labels}
  \item{col, lcol}{colour for plotting symbol and lines respectively}
  \item{colors}{vector of colours for each group}
  \item{alphalo, alphahi}{minimum and maximum transparency}
  \item{...}{other graphics parameters}
}
  

\value{
  Plot of density estimates (and partition) for discriminant analysis is
  sent to graphics window.
}

\details{
  If \code{prior.prob} is set to a particular value then this is used. 
  The default is \code{NULL} which means that the sample proportions are used.
  
  If \code{display="part"} then a partition induced by the discriminant
  analysis is also plotted. If this is not desired then set
  \code{display=""}. Its colours are controlled by \code{col}
  (the default is 2 to \eqn{\nu+1}{nu+1}, where \eqn{\nu}{nu} is the
  number of groups).   
  
  Unlike \code{\link{plot.kde}}, the contour plots are automatically
  added to the plot.  The line types are set by \code{lty} (the default
  is 1 to \eqn{\nu}{nu}). Also,
  \code{cont} and \code{ncont} control the number of level curves (only
  one of these needs to be set).

  The object \code{fhat} contains the training data and its group
  labels.  If \code{y} and \code{y.group} are missing then the training
  data points are plotted.  Otherwise, the test data \code{y} are plotted.
  The plotting symbols are set by \code{pch} (the default is 1 to
  \eqn{\nu}{nu}), one for each group.
}
 
\references{  

  Bowman, A.W. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Clarendon Press. Oxford.
  
  Simonoff, J. S., (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.}
}


\seealso{\code{\link{kda.kde}}, \code{\link{pda.pde}},
  \code{\link{kda}}, \code{\link{pda}}}

\examples{
## bivariate example

library(MASS)
data(iris)
ir <- iris[,c(1,2)]
ir.gr <- iris[,5]
xlab <- "Sepal length (mm)"
ylab <- "Sepal width (mm)"
xlim <- c(4,8)
ylim <- c(2,4.5)

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
fhat <- kda.kde(ir, ir.gr, H)
lda.fhat <- pda.pde(ir, ir.gr, type="line")
qda.fhat <- pda.pde(ir, ir.gr, type="quad")

layout(rbind(c(1,2), c(3,4)))
plot(fhat, cont=0, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, 
    pch=c(1,5,10))
plot(fhat, ncont=6, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
     col=c("transparent", "grey", "#8f8f8f"), drawlabels=FALSE)
plot(lda.fhat, ncont=6, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
     disp="")
plot(qda.fhat, ncont=6, xlim=xlim, ylim=ylim,  xlab=xlab, ylab=ylab, 
     lty=c(2,5,3))
layout(1)

## trivariate example
\dontrun{
ir <- iris[,1:3]
ir.gr <- iris[,5] 
H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
fhat <- kda.kde(ir, ir.gr, H)
plot(fhat, cont=c(25,50))
   ## colour indicates species, transparency indicates density heights

qda.fhat <- pda.pde(ir, ir.gr, type="quad")
plot(qda.fhat, cont=c(25,50))
}
}
\keyword{ hplot}

