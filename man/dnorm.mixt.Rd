\name{dnorm.mixt, rnorm.mixt}
\alias{rnorm.mixt}
\alias{dnorm.mixt}

\title{Univariate normal mixture distribution}
\description{
Random generation and density values from univariate normal mixture
distribution.}
\usage{
dnorm.mixt(x, mus=0, sigmas=1, props=1)
rnorm.mixt(n=100, mus=0, sigmas=1, props=1)
}

\arguments{
  \item{n}{number of random variates}
  \item{x}{vector of quantiles}
  \item{mus}{vector of means}
  \item{sigmas}{vector of standard deviations}
  \item{props}{vector of mixing proportions}
  }
  

\value{Univariate normal mixture random vectors and density values.}



\examples{
x <- rnorm.mixt(1000, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))
dens <- dnorm.mixt(x, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))
}

\keyword{ distribution }
