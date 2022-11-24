\name{vkde}
\alias{kde.balloon}
\alias{kde.sp}
\title{Variable kernel density estimate.}

\description{
  Variable kernel density estimate for 2-dimensional data.
}

\usage{
kde.balloon(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
   binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE)
kde.sp(x, H, h, gridsize, gridtype, xmin, xmax, supp=3.7, eval.points, 
   binned, bgridsize, w, compute.cont=TRUE, approx.cont=TRUE, verbose=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix. If this missing, \code{Hns}
  is called by default.}
  \item{h}{not yet implemented}
  \item{gridsize}{vector of number of grid points}
  \item{gridtype}{not yet implemented}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{supp}{effective support for standard normal}
  \item{eval.points}{vector or matrix of points at which estimate is evaluated}
  \item{binned}{flag for binned estimation.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{w}{vector of weights. Default is a vector of all ones.}
  \item{compute.cont}{flag for computing 1\% to 99\% probability contour levels. Default is TRUE.}
  \item{approx.cont}{flag for computing approximate probability contour
    levels. Default is TRUE.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
}
  
\value{
  A variable kernel density estimate for bounded data is an object of class \code{kde}.    
}

\details{
  The balloon density estimate \code{kde.balloon} employs bandwidths
  which vary at each 
  estimation point (Loftsgaarden & Quesenberry, 1965). There are as many bandwidths as there are estimation
  grid points. The default bandwidth is \code{Hns(,deriv.order=2)} and
  the subsequent bandwidths are derived via a minimal MSE formula.  

  The sample point density estimate \code{kde.sp} employs bandwidths
  which vary for each data point (Abramson, 1982).
  There are as many bandwidths as there are data
  points. The default bandwidth is \code{Hns(,deriv.order=4)} and the 
  subsequent bandwidths are derived via the Abramson formula.
}

\references{
  Abramson, I. S. (1982) On bandwidth variation in kernel estimates - a
  square root law. \emph{Annals of Statistics}, \bold{10}, 1217-1223.

  Loftsgaarden, D. O. & Quesenberry, C. P. (1965) A nonparametric
  estimate of a multivariate density function. \emph{Annals of
  Mathematical Statistics}, \bold{36}, 1049-1051.  
}

\seealso{\code{\link{kde}}, \code{\link{plot.kde}}}

\examples{
\donttest{data(worldbank)
wb <- as.matrix(na.omit(worldbank[,4:5]))
xmin <- c(-70,-35); xmax <- c(35,70)
fhat <- kde(x=wb, xmin=xmin, xmax=xmax)
fhat.sp <- kde.sp(x=wb, xmin=xmin, xmax=xmax)
zmax <- max(fhat.sp$estimate)
plot(fhat, display="persp", box=TRUE, phi=20, thin=1, border=grey(0,0.2), zlim=c(0,zmax))
plot(fhat.sp, display="persp", box=TRUE, phi=20, thin=1, border=grey(0,0.2), zlim=c(0,zmax))}
\dontrun{
fhat.ball <- kde.balloon(x=wb, xmin=xmin, xmax=xmax)
plot(fhat.ball, display="persp", box=TRUE, phi=20, zlim=c(0,zmax))}
}

\keyword{smooth}
