\name{pre.transform}
\alias{pre.sphere}
\alias{pre.scale}
\title{Pre-sphering and pre-scaling}

\description{
Pre-sphered or pre-scaled version of data.}

\usage{
pre.sphere(x, mean.centred=FALSE)
pre.scale(x, mean.centred=FALSE)
}

\arguments{
  \item{x}{matrix of data values}
  \item{mean.centred}{flag to centre the data values to have zero mean. Default is FALSE.}
}

\value{Pre-sphered or pre-scaled version of data. These
  pre-transformations are required for implementing the plug-in
  \code{\link{Hpi}} selectors and the smoothed cross validation
  \code{\link{Hscv}} selectors. 
}

\details{ For pre-scaling, the data values are pre-multiplied by
  \eqn{\mathbf{S}^{-1/2}}{S^(-1/2)} and for pre-scaling, by
  \eqn{\mathbf{S}_D^{-1/2}}{S_D^(-1/2)} where
  \eqn{\mathbf{S}}{S} is the sample variance and \eqn{\mathbf{S}_D}{S_D}
  is \eqn{\mathrm{diag} \, (S_1^2, S_2^2, \dots, S_d^2)}{diag (S_1^2,
    S_2^2, ..., S_d^2)} where
  \eqn{S_i^2}{S_i^2} is the i-th marginal sample variance. 
}

\examples{
data(unicef)
unicef.sp <- pre.sphere(as.matrix(unicef))
}

\keyword{algebra}
