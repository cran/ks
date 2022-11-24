\name{kdcde}
\alias{kdcde}
\alias{dckde}
\title{Deconvolution kernel density derivative estimate}

\description{
  Deconvolution kernel density derivative estimate for 1- to 6-dimensional data.
}

\usage{
kdcde(x, H, h, Sigma, sigma, reg, bgridsize, gridsize, binned, 
      verbose=FALSE, ...)
dckde(...) 
}

\arguments{
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth matrix/scalar bandwidth. If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{Sigma,sigma}{error variance matrix}
  \item{reg}{regularisation parameter}
  \item{gridsize}{vector of number of grid points}
  \item{binned}{flag for binned estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is
    FALSE.}
  \item{...}{other parameters to \code{\link{kde}}} 
}
  
\value{
  A deconvolution kernel density derivative estimate is an object of class
  \code{kde} which is a list with fields:
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{h}{scalar bandwidth (1-d only)}
  \item{H}{bandwidth matrix}
  \item{gridtype}{"linear"}
  \item{gridded}{flag for estimation on a grid}
  \item{binned}{flag for binned estimation}
  \item{names}{variable names} 
  \item{w}{vector of weights}
  \item{cont}{vector of probability contour levels}   
}

\details{ 
  A weighted kernel density estimate is utilised to perform the
  deconvolution. The weights \code{w} are the solution to a
  quadratic programming problem, and then input into \code{kde(,w=w)}.
  This weighted estimate also requires an estimate of the error
  variance matrix from repeated observations, and of the regularisation
  parameter. If the latter is missing, it is calculated internally using
  a 5-fold cross validation method. See Hazelton & Turlach (2009).
  \code{dckde} is an alias for \code{kdcde}.
  
  If the bandwidth \code{H} is missing from \code{kde}, then
  the default bandwidth is the plug-in selector
  \code{Hpi}. Likewise for missing \code{h}.
  
  The effective support, binning, grid size, grid range, positive
  parameters are the same as \code{\link{kde}}.
}

\references{
Hazelton, M. L. & Turlach, B. A. (2009), Nonparametric density
deconvolution by weighted kernel density estimators, \emph{Statistics
  and Computing}, \bold{19}, 217-228.
}

\seealso{\code{\link{kde}}}

\examples{
\donttest{
data(air)
air <- air[, c("date", "time", "co2", "pm10")]
air2 <- reshape(air, idvar="date", timevar="time", direction="wide")
air <- as.matrix(na.omit(air2[,c("co2.20:00", "pm10.20:00")]))
Sigma.air <- diag(c(var(air2[,"co2.19:00"] - air2["co2.21:00"], na.rm=TRUE),
   var(air2[,"pm10.19:00"] - air2[,"pm10.21:00"], na.rm=TRUE)))
fhat.air.dec <- kdcde(x=air, Sigma=Sigma.air, reg=0.00021, verbose=TRUE)
plot(fhat.air.dec, drawlabels=FALSE, display="filled.contour", lwd=1)}
}

\keyword{smooth}
