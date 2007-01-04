\name{dmvnorm.mixt, rmvnorm.mixt}
\alias{rmvnorm.mixt}
\alias{dmvnorm.mixt}


\title{Multivariate normal mixture distribution}
\description{
Random generation and density values from multivariate normal mixture
distribution.}
\usage{
rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1)
dmvnorm.mixt(x, mus, Sigmas, props)
}

\arguments{
  \item{n}{number of random variates}
  \item{x}{matrix of quantiles}
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  }
  

\value{Multivariate  normal mixture random vectors and density values.}

\details{
  \code{rmvnorm.mixt} and \code{dmvnorm.mixt} are based on the
  \code{rmvnorm} and \code{dmvnorm} functions from the \code{mvtnorm} library.
}  

\seealso{\code{\link{rmvt.mixt}}, \code{\link{dmvt.mixt}},
  \code{\link{dnorm.mixt}}, \code{\link{rnorm.mixt}}}

\examples{
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(1000, mus, Sigmas, props)
dens <- dmvnorm.mixt(x, mus, Sigmas, props)
}

\keyword{ distribution }
