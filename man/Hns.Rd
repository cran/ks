\name{Hns}
\alias{Hns}
\alias{hns}

\title{Normal scale bandwidth}

\description{
Normal scale bandwidth.  
}


\usage{
Hns(x, deriv.order=0)
hns(x, deriv.order=0)
}

\arguments{
  \item{x}{vector/matrix of data values}
  \item{deriv.order}{derivative order}
}


\value{
  Full normal scale bandwidth matrix.  
}

\details{
  \code{Hns} is equal to \code{(4/(n*(d+2*r+2)))^(2/(d+2*r+4))*var(x)}, n = sample size, d = dimension of data, r = derivative order. \code{hns} is the analogue  of \code{Hns} for 1-d data.  
}


\references{Chacon J.E., Duong, T. & Wand, M.P. (2011). Asymptotics for
	general multivariate kernel density derivative
	estimators. \emph{Statistica Sinica}. \bold{21}, 807-840.  
}

\examples{
x <- rmvnorm.mixt(1000)
Hns(x, deriv.order=1)
hns(x[,1], deriv.order=1)
}

\keyword{smooth}

