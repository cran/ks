\name{Hns}
\alias{Hns}
\alias{Hns.diag}
\alias{hns}
\alias{Hns.kcde}
\alias{hns.kcde}
\title{Normal scale bandwidth}

\description{
Normal scale bandwidth.  
}

\usage{
Hns(x, deriv.order=0)
Hns.diag(x)
hns(x, deriv.order=0)
Hns.kcde(x)
hns.kcde(x)
}

\arguments{
  \item{x}{vector/matrix of data values}
  \item{deriv.order}{derivative order}
}

\value{
  Normal scale bandwidth. 
}

\details{
  \code{Hns} is equal to \code{(4/(n*(d+2*r+2)))^(2/(d+2*r+4))*var(x)},
  n = sample size, d = dimension of data, r = derivative
  order. \code{hns} is the analogue  of \code{Hns} for 1-d data. These
  can be used for density (derivative) estimators
  \code{\link{kde}}, \code{\link{kdde}}. 
  The equivalents for distribution estimators \code{\link{kcde}} are
  \code{Hns.kcde} and \code{hns.cde}. 
}

\references{Chacon J.E., Duong, T. & Wand, M.P. (2011). Asymptotics for
  general multivariate kernel density derivative
  estimators. \emph{Statistica Sinica}, \bold{21}, 807-840.   
}

\examples{
data(forbes, package="MASS")
Hns(forbes, deriv.order=2)
hns(forbes$bp, deriv.order=2)
}

\keyword{smooth}

