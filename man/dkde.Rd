\name{dkde}
\alias{dkde}
\alias{pkde}
\alias{qkde}
\alias{rkde}

\title{Functions for univariate kernel density estimates}

\description{
   Functions for 1-dimensional kernel density estimates.
}


\usage{
 pkde(q, fhat)
 qkde(p, fhat)
 dkde(x, fhat)
 rkde(n, fhat, positive=FALSE) 
}

\arguments{
   \item{x,q}{vector of quantiles}
   \item{p}{vector of probabilities}
   \item{n}{number of observations}
   \item{positive}{flag to compute KDE on the positive real
	 line. Default is FALSE.} 
   \item{fhat}{kernel density estimate, object of class \code{"kde"}}
}


\value{
   For the kernel density estimate \code{fhat},
   \code{pkde} computes the cumulative probability for the quantile
   \code{q}, \code{qkde} computes the quantile corresponding to the probability
   \code{p}, \code{dkde} computes the density value at
   \code{x} and \code{rkde} computes a random sample of size \code{n}.  
}

\details{ 
  \code{pkde} uses Simpson's rule for the numerical
  integration. \code{rkde} uses
  Silverman (1986)'s method to generate a random sample from a KDE.
}

\references{ Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman & Hall/CRC. London.}

\examples{
x <- rnorm.mixt(n=10000, mus=0, sigmas=1, props=1)
fhat <- kde(x=x, h=hpi(x))
p1 <- pkde(fhat=fhat, q=c(-1, 0, 0.5))
qkde(fhat=fhat, p=p1)    
y <- rkde(fhat=fhat, n=100)

## positive data example
z <- 2^x
fhatz <- kde(x=z, h=hpi(log(z)), positive=TRUE)
z2 <- rkde(fhat=fhatz, n=10000, positive=TRUE)
fhatz2 <- kde(x=z2, h=hpi(log(z2)), positive=TRUE)
plot(fhatz, xlim=c(0,10))
plot(fhatz2, add=TRUE, col=2) 
}

\keyword{ smooth}

