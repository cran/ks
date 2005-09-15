\name{kde}
\alias{kde}

\title{Kernel density estimate for multivariate data}
\description{
  Kernel density estimate for 2- to 6-dimensional data
}
\usage{
 kde(x, H, gridsize, supp=3.7, eval.points)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{gridsize}{vector of number of grid points}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points that density estimate is evaluated at (required
       for dimensions > 3)} 
}
  

\value{
  Kernel density estimate is an object of class \code{kde} which is a
  list with 4 fields
  \item{x}{data points - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrix}
}

\details{The kernel density estimate is computed exactly i.e. binning
  is not used. 
  
  If \code{gridsize} is not set to a specific value, then it
  defaults to 50 grid points in each co-ordinate direction
  i.e. \code{rep(50, d)}. Not required
  to be set if specifying \code{eval.points}.

  If \code{eval.points} is not specified, then the
  density estimate is automatically computed over a grid whose
  resolution is controlled by \code{gridsize} (a grid is
  required for plotting). 
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\seealso{\code{\link{plot.kde}}}

\examples{
### bivariate example
data(unicef)
H.pi <- Hpi(unicef, nstage=1)
fhat <- kde(unicef, H.pi)


### 4-variate example
library(MASS)
data(iris)
ir <- iris[,1:4][iris[,5]=="setosa",]
H.scv <- Hscv(ir)
fhat <- kde(ir, H.scv, eval.points=ir)  
}

\keyword{smooth}
