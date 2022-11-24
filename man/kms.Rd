\name{kms}
\alias{kms}
\alias{summary.kms}
\alias{plot.kms}
\title{Kernel mean shift clustering}

\description{
  Kernel mean shift clustering for 2- to 6-dimensional data.
}

\usage{
kms(x, y, H, max.iter=400, tol.iter, tol.clust, min.clust.size, merge=TRUE,
    keep.path=FALSE, verbose=FALSE)

\method{plot}{kms}(x, display="splom", col, col.fun, alpha=1, xlab, ylab, zlab, theta=-30, 
    phi=40, add=FALSE, ...)
\method{summary}{kms}(object, ...)
}

\arguments{
  \item{x}{matrix of data values or object of class \code{kms}}
  \item{y}{matrix of candidate data values for which the mean shift will
    estimate their cluster labels. If missing, \code{y=x}.}
  \item{H}{bandwidth matrix/scalar bandwidth. If missing,
  \code{Hpi(x,deriv.order=1,nstage=2-(d>2))} is called by default.}
  \item{max.iter}{maximum number of iterations. Default is 400.}
  \item{tol.iter}{distance under which two successive iterations are
    considered convergent. Default is 0.001*min marginal IQR of \code{x}.}
  \item{tol.clust}{distance under which two cluster modes are considered
    to form one cluster. Default is 0.01*max marginal IQR of \code{x}.}
  \item{min.clust.size}{minimum cluster size (cardinality). Default is \code{0.01*nrow(y)}.}
  \item{merge}{flag to merge clusters which are smaller than
    \code{min.clust.size}. Default is TRUE.}
  \item{keep.path}{flag to store the density gradient ascent paths. Default is FALSE.}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{object}{object of class \code{kms}}
  \item{display}{type of display, "splom" (>=2-d) or "plot3D" (3-d)}
  \item{col,col.fun}{vector or colours (one for each group) or colour function}
  \item{alpha}{colour transparency. Default is 1.}
  \item{xlab,ylab,zlab}{axes labels}
  \item{theta,phi}{graphics parameters for perspective plots (3-d)}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{...}{other (graphics) parameters}
}
  
\value{
  A kernel mean shift clusters set is an object of class \code{kms}
  which is a list with fields: 
  \item{x,y}{data points - same as input}
  \item{end.points}{matrix of final iterates starting from \code{y}}
  \item{H}{bandwidth matrix}
  \item{label}{vector of cluster labels}
  \item{nclust}{number of clusters}
  \item{nclust.table}{frequency table of cluster labels}
  \item{mode}{matrix of cluster modes}
  \item{names}{variable names}
  \item{tol.iter,tol.clust,min.clust.size}{tuning parameter values -
  same as input}
  \item{path}{list of density gradient ascent paths where \code{path[[i]]} is
    the path of \code{y[i,]} (only if \code{keep.path=TRUE})}
}

\details{ 
  Mean shift clustering belongs to the class of modal or density-based
  clustering methods.
  The mean shift recurrence of the candidate point \eqn{{\bold x}}{x} is
  \eqn{{\bold x}_{j+1} = {\bold x}_j + \bold{{\rm H}} {\sf D} \hat{f}({\bold
      x}_j)/\hat{f}({\bold x}_j)}{x_j+1 = x_j + H D hat(f)(x_j)/hat(f)(x_j)}
  where \eqn{j\geq 0}{j>=0} and \eqn{{\bold x}_0 = {\bold x}}{x_0 = x}.
  The sequence \eqn{\{{\bold x}_0, {\bold
  x}_1, \dots \}}{x_0, x_1, ...} follows the density gradient ascent
  paths to converge to a local mode of the
  density estimate \eqn{\hat{f}}{hat(f)}. Hence \eqn{{\bold x}}{x} is
  iterated until it converges to its local mode, and this determines its
  cluster label.  

  The mean shift recurrence is terminated if successive iterations are
  less than \code{tol.iter} or the maximum number of iterations
  \code{max.iter} is reached. Final iterates which are less than
  \code{tol.clust} distance apart are considered to form a single
  cluster. If \code{merge=TRUE} then the clusters whose cardinality is less
  than \code{min.clust.size} are iteratively merged with their nearest cluster.   
  
  If the bandwidth \code{H} is missing, then
  the default bandwidth is the plug-in selector for the density gradient
  \code{Hpi(x,deriv.order=1)}. Any bandwidth that is suitable for the
  density gradient is also suitable for the mean shift. 
}

\references{
  Chacon, J.E. & Duong, T. (2013) Data-driven density estimation, with
  applications to nonparametric clustering and bump hunting. \emph{Electronic
    Journal of Statistics}, \bold{7}, 499-532.

  Comaniciu, D. & Meer, P. (2002). Mean shift: a robust approach
  toward feature space analysis. \emph{ IEEE Transactions on Pattern
  Analysis and Machine Intelligence}, \bold{24}, 603-619.
}  

\seealso{\code{\link{kde}}}

\examples{
library(MASS)
data(crabs)
kms.crabs <- kms(x=crabs[,c("FL","CW")])
plot(kms.crabs, pch=16)
summary(kms.crabs)

kms.crabs <- kms(x=crabs[,c("FL","CW","RW")])
plot(kms.crabs, pch=16)
plot(kms.crabs, display="plot3D", pch=16) 
}

\keyword{cluster}
