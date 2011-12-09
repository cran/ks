\name{dmvnorm.mixt}
\alias{rnorm.mixt}
\alias{dnorm.mixt}
\alias{rmvnorm.mixt}
\alias{dmvnorm.mixt}


\title{Multivariate and univariate normal mixture distribution}
\description{
Random generation and density values from normal mixture distribution.}
\usage{
dmvnorm.mixt(x, mus, Sigmas, props)
rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1, mixt.label=FALSE)

dnorm.mixt(x, mus=0, sigmas=1, props=1)
rnorm.mixt(n=100, mus=0, sigmas=1, props=1, mixt.label=FALSE)
}

\arguments{
  \item{n}{number of random variates}
  \item{x}{matrix of quantiles}
  \item{mus}{(stacked) matrix of mean vectors (>1-d) or vector of means (1-d)}
  \item{Sigmas}{(stacked) matrix of variance matrices (>1-d)}
  \item{sigmas}{vector of standard deviations (1-d)}
  \item{props}{vector of mixing proportions}
  \item{mixt.label}{flag to output numeric label indicating mixture component. Default is FALSE.}
  }
  

\value{Normal mixture random vectors and density values.}

\details{
  \code{rmvnorm.mixt} and \code{dmvnorm.mixt} are based on the
  \code{rmvnorm} and \code{dmvnorm} functions from the \code{mvtnorm} library.
}  

\seealso{\code{\link{rmvt.mixt}}, \code{\link{dmvt.mixt}},
  \code{\link{dnorm.mixt}}, \code{\link{rnorm.mixt}}}

\examples{
## bivariate 
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(1000, mus, Sigmas, props)
dens <- dmvnorm.mixt(x, mus, Sigmas, props)

##univariate
x <- rnorm.mixt(1000, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))
dens <- dnorm.mixt(x, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))
}

\keyword{ distribution }
