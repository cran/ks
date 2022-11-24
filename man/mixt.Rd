\name{mixt}
\alias{rnorm.mixt}
\alias{dnorm.mixt}
\alias{rmvnorm.mixt}
\alias{dmvnorm.mixt}
\alias{rmvt.mixt}
\alias{dmvt.mixt}
\alias{mvnorm.mixt.mode}
\title{Normal and t-mixture distributions}

\description{
Random generation and density values from normal and t-mixture distributions.}

\usage{
dnorm.mixt(x, mus=0, sigmas=1, props=1)
rnorm.mixt(n=100, mus=0, sigmas=1, props=1, mixt.label=FALSE)
dmvnorm.mixt(x, mus, Sigmas, props=1, verbose=FALSE)
rmvnorm.mixt(n=100, mus=c(0,0), Sigmas=diag(2), props=1, mixt.label=FALSE)
rmvt.mixt(n=100, mus=c(0,0), Sigmas=diag(2), dfs=7, props=1)
dmvt.mixt(x, mus, Sigmas, dfs, props)
mvnorm.mixt.mode(mus, Sigmas, props=1, verbose=FALSE)
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
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{dfs}{vector of degrees of freedom}
}
  
\value{Normal and t-mixture random vectors and density values.}

\details{
  \code{rmvnorm.mixt} and \code{dmvnorm.mixt} are based on the
  \code{rmvnorm} and \code{dmvnorm} functions from the \pkg{mvtnorm}
  package. Likewise for \code{rmvt.mixt} and \code{dmvt.mixt}.

  For the normal mixture densities, \code{mvnorm.mixt.mode} computes the
  local modes: these are usually very close but not exactly equal to the
  component means.
}

\examples{
## univariate normal mixture
x <- rnorm.mixt(1000, mus=c(-1,1), sigmas=c(0.5, 0.5), props=c(1/2, 1/2))

## bivariate mixtures 
mus <- rbind(c(-1,0), c(1, 2/sqrt(3)), c(1,-2/sqrt(3)))
Sigmas <- 1/25*rbind(invvech(c(9, 63/10, 49/4)), invvech(c(9,0,49/4)), invvech(c(9,0,49/4)))
props <- c(3,3,1)/7
dfs <- c(7,3,2)
x <- rmvnorm.mixt(1000, mus=mus, Sigmas=Sigmas, props=props)
y <- rmvt.mixt(1000, mus=mus, Sigmas=Sigmas, dfs=dfs, props=props)

mvnorm.mixt.mode(mus=mus, Sigmas=Sigmas, props=props)
}

\keyword{distribution}
