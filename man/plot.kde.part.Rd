\name{plot.kde.part}
\alias{plot.kde.part}
\alias{kms.part}
\alias{mvnorm.mixt.part}
\title{Partition plot for kernel density clustering}

\description{
  Plot of partition for kernel density clustering for 2-dimensional data.
}

\usage{
mvnorm.mixt.part(mus, Sigmas, props=1, xmin, xmax, gridsize, max.iter=100,
   verbose=FALSE)
kms.part(x, H, xmin, xmax, gridsize, verbose=FALSE, ...)

\method{plot}{kde.part}(x, display="filled.contour", col, col.fun, alpha=1, add=FALSE, ...)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{gridsize}{vector of number of grid points}
   \item{max.iter}{maximum number of iterations}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{x}{matrix of data values or an object of class \code{kde.part}}
  \item{H}{bandwidth matrix. If missing,
    \code{Hpi(x,deriv,order=1)} is called by default.}
  \item{display}{type of display, "filled.contour" for filled contour plot} 
  \item{col,col.fun}{vector of plotting colours or colour function}
  \item{alpha}{colour transparency. Default is 1.}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{...}{other parameters}
}
  
\value{
  A kernel partition is an object of class \code{kde.part} which is a
  list with fields:
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{w}{vector of weights}
  \item{cont}{vector of probability contour levels}   
  \item{end.points}{matrix of final iterates starting from \code{x}}
  \item{label}{vector of cluster labels}
  \item{mode}{matrix of cluster modes}
  \item{nclust}{number of clusters}
  \item{nclust.table}{frequency table of cluster labels}
  \item{tol.iter,tol.clust,min.clust.size}{tuning parameter values -
    same as input}
  
  Plot is sent to graphics window. 
}

\details{
  For 2-d data, \code{kms.part} and \code{mvnorm.mixt.part} produce a
  \code{kde.part} object whose 
  values are the class labels, rather than probability density values.   
}

\seealso{\code{\link{plot.kde}}, \code{\link{kms}}}

\examples{
## normal mixture partition
mus <- rbind(c(-1,0), c(1, 2/sqrt(3)), c(1,-2/sqrt(3)))
Sigmas <- 1/25*rbind(invvech(c(9, 63/10, 49/4)), invvech(c(9,0,49/4)), invvech(c(9,0,49/4)))
props <- c(3,3,1)/7
gridsize <- c(11,11) ## small gridsize illustrative purposes only 
nmixt.part <- mvnorm.mixt.part(mus=mus, Sigmas=Sigmas, props=props, gridsize=gridsize)
plot(nmixt.part, asp=1, xlim=c(-3,3), ylim=c(-3,3), alpha=0.5)

## kernel mean shift partition
\donttest{set.seed(81928192)
x <- rmvnorm.mixt(n=10000, mus=mus, Sigmas=Sigmas, props=props)
msize <- round(prod(gridsize)*0.1)
kms.nmixt.part <- kms.part(x=x, min.clust.size=msize, gridsize=gridsize)
plot(kms.nmixt.part, asp=1, xlim=c(-3,3), ylim=c(-3,3), alpha=0.5)}}

\keyword{hplot}
