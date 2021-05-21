\name{kdr}
\alias{kdr}
\alias{plot.kdr}

\title{Kernel density ridge estimation}
\description{
  Kernel density ridge estimation for 2- to 3-dimensional data.
}
\usage{
kdr(x, y, H, p=1, max.iter=400, tol.iter, tol.seg, min.seg.size,
    keep.path=FALSE, gridsize, xmin, xmax, binned, bgridsize, w, fhat,
    density.cutoff, pre=TRUE, verbose=FALSE) 
    
\method{plot}{kdr}(x, add=FALSE, ...)
}

\arguments{
  \item{x}{matrix of data values or an object of class \code{kdr}}
  \item{y}{matrix of initial values}
  \item{p}{dimension of density ridge}
  \item{H}{bandwidth matrix/scalar bandwidth. If missing,
  \code{Hpi(x,deriv,order=2)} is called by default.}
  \item{max.iter}{maximum number of iterations. Default is 400.}
  \item{tol.iter}{distance under which two successive iterations are
    considered convergent. Default is 0.001*min marginal IQR of \code{x}.}
  \item{tol.seg}{distance under which two segments are considered
    to form one segment. Default is 0.01*max marginal IQR of \code{x}.}
  \item{min.seg.size}{minimum length of a segment of a density
  ridge. Default is \code{round(0.001*nrow(y),0)}.}
  \item{keep.path}{flag to store the density gradient ascent
    paths. Default is FALSE.}
  \item{gridsize}{vector of number of grid points}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{binned}{flag for binned estimation.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{fhat}{kde of \code{x}. If missing \code{kde(x=x,w=w)} is
    executed.}
  \item{density.cutoff}{density threshold under which the \code{y} are
    excluded from the density ridge estimation. Default is
  \code{contourLevels(fhat, cont=99)}.}
  \item{pre}{flag for pre-scaling data. Default is TRUE.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{...}{other graphics parameters}
}
  

\value{
  A kernel density ridge set is an object of class \code{kdr}
  which is a list with fields: 
  \item{x,y}{data points - same as input}
  \item{end.points}{matrix of final iterates starting from \code{y}}
  \item{H}{bandwidth matrix}
  \item{names}{variable names}
  \item{tol.iter,tol.clust,min.seg.size}{tuning parameter values -
    same as input}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{w}{vector of weights}
  \item{path}{list of density gradient ascent paths where \code{path[[i]]} is
    the path of \code{y[i,]} (only if \code{keep.path=TRUE})}
}

\details{ 
  Kernel density ridge estimation is based on reduced dimension kernel
  mean shift. See Ozertem & Erdogmus (2011). 

  If \code{y} is missing, then it defaults to the grid of size
  \code{gridsize} spanning from \code{xmin} to \code{xmax}.
  
  If the bandwidth \code{H} is missing, then
  the default bandwidth is the plug-in selector for the density gradient
  \code{Hpi(x,deriv.order=2)}. Any bandwidth that is suitable for the
  density Hessian is also suitable for the kernel density ridge. 
}

\references{
  Ozertem, U. & Erdogmus, D. (2011), Locally defined principal curves and
  surfaces, \emph{Journal of Machine Learning Research}, \bold{12}, 
  1249-1286. 
}  

\seealso{\code{\link{kde}}}

\examples{
\donttest{library(maps)
data(quake) 
quake <- quake[quake$prof==1,]  ## Pacific Ring of Fire 
quake$long[quake$long<0] <- quake$long[quake$long<0] + 360
quake <- quake[, c("long", "lat")]
data(plate)                     ## tectonic plate boundaries
plate <- plate[plate$long < -20 | plate$long > 20,]
plate$long[plate$long<0 & !is.na(plate$long)] <- plate$long[plate$long<0
& !is.na(plate$long)] + 360

dr.quake <- kdr(x=quake, xmin=c(70,-70), xmax=c(310, 80))
map(wrap=c(0,360), lty=2)
lines(plate[,1:2], col=4, lwd=2)
plot(dr.quake, cex=0.5, pch=16, col=6, add=TRUE)}
}

\keyword{cluster}
