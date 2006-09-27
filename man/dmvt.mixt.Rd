\name{dmvt.mixt, rmvt.mixt}
\alias{rmvt.mixt}
\alias{dmvt.mixt}

\title{Multivariate t mixture distribution}
\description{
Random generation and density values from multivariate t mixture 
distribution.}
\usage{
rmvt.mixt(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
dmvt.mixt(x, mus, Sigmas, dfs, props)
}

\arguments{
  \item{n}{number of random variates}
  \item{x}{matrix of quantiles}
  \item{mus}{(stacked) matrix of location vectors}
  \item{Sigmas}{(stacked) matrix of dispersion matrices}
  \item{dfs}{vector of degrees of freedom}
  \item{props}{vector of mixing proportions}
  }
  

\value{Multivariate t mixture random vectors and density values.}

\details{
  \code{rmvt.mixt} and \code{dmvt.mixt} are based on the
  \code{rmvt} and \code{dmvt} functions from the \code{mvtnorm} library.

  The formula for a d-variate t density with location vector
  \eqn{\bold{\mu}}{mu}, dispersion matrix
  \eqn{\bold{\Sigma}}{Sigma} and df degrees of freedom is
  \deqn{\frac{\Gamma((df+d)/2)}{(df\pi)^{d/2} \Gamma(df/2)
      |\bold{\Sigma}^{1/2}|}\left[ 1 +  \frac{1}{df} (\bold{x} - \bold{\mu})^T
    \bold{\Sigma}^{-1}(\bold{x} - \bold{\mu})
    \right]^{-(d+df)/2}.}{gamma((df+d)/2) / ((df*pi)^(d/2) * gamma(df/2)* det(Sigma)^(1/2) * (1 + 1/df * (x-mu)^T * Sigma^(-1) (x-mu))^(-(d+df)/2).} 
}


\seealso{\code{\link{rmvnorm.mixt}}, \code{\link{dmvnorm.mixt}}}

\examples{
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
dfs <- c(7,3)
x <- rmvt.mixt(1000, mus, Sigmas, dfs, props)
dens <- dmvt.mixt(x, mus, Sigmas, dfs, props)
}

\keyword{ distribution  }

