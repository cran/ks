\name{pre.scale, pre.sphere}
\alias{pre.sphere}
\alias{pre.scale}

\title{Pre-sphering and pre-scaling}
\description{
Pre-sphered or pre-scaled version of data.}
\usage{
pre.sphere(x)
pre.scale(x)
}

\arguments{
  \item{x}{matrix of data values}
}
\value{Pre-sphered or pre-scaled version of data. These
  pre-transformations are required for implementing the plug-in
  \code{\link{Hpi}} selectors and the smoothed cross validation
  \code{\link{Hscv}} selectors. 
}

\details{ For pre-scaling, the data values are pre-multiplied by
  \eqn{\mathbf{S}^{-1/2}}{S^(-1/2)} and for pre-scaling, by
  \eqn{(\mathbf{S}_D)^{-1/2}}{S_D^(-1/2)} where
  \eqn{\mathbf{S}}{S} is the sample variance and \eqn{\mathbf{S}_D}{S_D}
  is \eqn{\mathrm{diag} \, (S_1^2, S_2^2, \dots, S_d^2)}{diag (S_1^2,
    S_2^2, ..., S_d^2)} where
  \eqn{S_i^2}{S_i^2} is the i-th marginal sample variance.

  If \eqn{\mathbf{H}^*}{H*} is the bandwidth matrix for the
  pre-transformed data and \eqn{\mathbf{H}}{H} is the bandwidth matrix for the
  original data, then
  \eqn{\mathbf{H}=\mathbf{S}^{1/2} \mathbf{H}^* \mathbf{S}^{1/2}}{S^(1/2) H* S^(1/2)} or \eqn{\mathbf{H} = \mathbf{S}_D^{1/2} \mathbf{H}^*\mathbf{S}_D^{1/2}}{S_D^(1/2) H* S_D^(1/2)} as appropriate.
  
}

\references{ Wand, M.P. \& Jones, M.C. (1994) Multivariate plugin bandwidth
    selection. \emph{Computational Statistics}, \bold{9}, 97-116.
  
  Duong, T. \& Hazelton, M.L. (2003) Plug-in bandwidth matrices for
    bivariate kernel density estimation. \emph{Journal of Nonparametric
  Statistics}, \bold{15}, 17-30.
}

\examples{
x <- rmvnorm.mixt(1000, mus=c(0,0), Sigmas=rbind(c(1,0.2), c(0.2, 0.5))) 
x.sp <- pre.sphere(x)
x.sc <- pre.scale(x)
var(x.sp)
var(x.sc)
}
\keyword{ algebra }
