\name{plot.kde}
\alias{plot.kde}

\title{Kernel density estimate plot for 2- and 3-dimensional data}
\description{
  Kernel density estimate plot for 2- and 3-dimensional data.
}

\usage{
\method{plot}{kde}(x, display="slice", ...)
}

\arguments{
  \item{x}{an object of class \code{kde} i.e. output from
    \code{\link{kde}} function}
  \item{display}{type of display}
  \item{...}{other graphics parameters - see details below}
}
  

\value{
  Plot of kernel density estimate is sent to graphics window.
}

\references{
  Bowman, A.W. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Clarendon Press. Oxford.
  
  Simonoff, J. S., (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.}
}

\details{  
  There are three types of plotting displays available, controlled
  by the \code{display} parameter.

  If \code{display="slice"} then a slice/contour plot
  is generated using \code{contour}.
  The default contours are at 25\%, 50\%, 75\% or
  \code{cont=c(25,50,75)}.  The user can also set the number of contour
  level curves by changing the value set to \code{ncont}. See examples below.
  
  If \code{display="persp"} then a perspective/wire-frame plot
  is generated.  The default z-axis limits \code{zlim} are determined by
  the range of the z values i.e. default from the usual \code{persp} command.
  
  If \code{display="image"} then an image plot
  is generated. The colors are the default from the usual
  \code{image} command.

  For 3-dimensional data, the plot is a series of 2-dimensional contour plots. 
  Use \code{layout.mat} to control the grid layout from the usual \code{layout}
  command.     
}
 

\seealso{\code{\link{kde}}}

\examples{
### bivariate example
data(unicef)
H.scv <- Hscv(unicef)
fhat <- kde(unicef, H.scv)

layout(rbind(c(1,2), c(3,4)))
plot(fhat, display="slice", cont=seq(10,90, by=20), cex=0.3)
plot(fhat, display="slice", ncont=5, cex=0.3, drawlabels=FALSE)
plot(fhat, display="persp")
plot(fhat, display="image", col=rev(heat.colors(15)))
layout(1)

### 3-variate example

mus <- rbind(c(0,0,0), c(2,2,2))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
x <- rmvnorm.mixt(n=100, mus=mus, Sigmas=Sigmas, props=props)
H.pi <- Hpi(x)
fhat <- kde(x, H.pi, eval.levels=seq(-3,3, length=9))  
plot(fhat, disp="slice", ncont=6, cex=0.3)
}

\keyword{ hplot}
