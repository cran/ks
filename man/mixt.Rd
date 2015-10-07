\name{mixt}
\alias{rnorm.mixt}
\alias{dnorm.mixt}
\alias{rmvnorm.mixt}
\alias{dmvnorm.mixt}
\alias{rmvt.mixt}
\alias{dmvt.mixt}

\title{Normal and t-mixture distributions}
\description{
Random generation and density values from normal and t-mixture distributions.}
\usage{
dnorm.mixt(x, mus=0, sigmas=1, props=1)
rnorm.mixt(n=100, mus=0, sigmas=1, props=1, mixt.label=FALSE)
dmvnorm.mixt(x, mus, Sigmas, props)
rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1, mixt.label=FALSE)
rmvt.mixt(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
dmvt.mixt(x, mus, Sigmas, dfs, props)
}

\arguments{
  \item{n}{number of random variates}
  \item{x}{matrix of quantiles}
  \item{mus}{(stacked) matrix of mean vectors (>1-d) or vector of means (1-d)}
  \item{Sigmas}{(stacked) matrix of variance matrices (>1-d)}
  \item{sigmas}{vector of standard deviations (1-d)}
  \item{props}{vector of mixing proportions}
  \item{mixt.label}{flag to output numeric label indicating mixture
    component. Default is FALSE.}
  \item{dfs}{vector of degrees of freedom}
  }
  

\value{Normal and t-mixture random vectors and density values.}

\details{
  \code{rmvnorm.mixt} and \code{dmvnorm.mixt} are based on the
  \code{rmvnorm} and \code{dmvnorm} functions from the \pkg{mvtnorm}
  package. Likewise for \code{rmvt.mixt} and \code{dmvt.mixt}.
}  

\examples{
## univariate normal mixture
x <- rnorm.mixt(1000, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))

## bivariate mixtures 
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
dfs <- c(7,3)
x <- rmvnorm.mixt(1000, mus=mus, Sigmas=Sigmas, props=props)
y <- rmvt.mixt(1000, mus=mus, Sigmas=Sigmas, dfs=dfs, props=props)
}

\keyword{ distribution }
