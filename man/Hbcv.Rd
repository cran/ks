\name{Hbcv}
\alias{Hbcv}
\alias{Hbcv.diag}
\title{Biased cross-validation (BCV) bandwidth matrix selector
  for bivariate data}

\description{
BCV bandwidth matrix for bivariate data.}

\usage{
Hbcv(x, whichbcv=1, Hstart, binned=FALSE, amise=FALSE, verbose=FALSE)
Hbcv.diag(x, whichbcv=1, Hstart, binned=FALSE, amise=FALSE, verbose=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{whichbcv}{1 = BCV1, 2 = BCV2.  See details below.}
  \item{Hstart}{initial bandwidth matrix, used in numerical optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{amise}{flag to return the minimal BCV value. Default is FALSE.}
  \item{verbose}{flag to print out progress information. Default is FALSE.}  
}

\value{
BCV bandwidth matrix. If \code{amise=TRUE} then the minimal BCV value is returned too. 
}

\references{Sain, S.R, Baggerly, K.A. & Scott, D.W. (1994)
  Cross-validation of multivariate densities. \emph{Journal of the
  American Statistical Association}, \bold{82}, 1131-1146.
}

\details{
  Use \code{Hbcv} for unconstrained bandwidth matrices and \code{Hbcv.diag}
  for diagonal bandwidth matrices. These selectors are only
  available for bivariate data. Two types of BCV criteria are
  considered here.  They are known as BCV1 and BCV2, from Sain, Baggerly
  & Scott (1994) and only differ slightly. These BCV
  surfaces can have multiple minima and so it can be quite difficult to
  locate the most appropriate minimum. Some times, there can be no local minimum at all so there
  may be no finite BCV selector.

  For details about the advanced options for \code{binned}, \code{Hstart}, see \code{\link{Hpi}}.
}

\seealso{\code{\link{Hlscv}}, \code{\link{Hpi}}, \code{\link{Hscv}}}

\examples{
data(unicef)
Hbcv(unicef)
Hbcv.diag(unicef)
}

\keyword{smooth}
