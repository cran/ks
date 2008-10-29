\name{dkde, pkde, qkde, rkde}
\alias{dkde}
\alias{pkde}
\alias{qkde}
 \alias{rkde}

 \title{Functions for 1-dimensional kernel density estimates}

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
  \code{pkde} uses the Simpson's rule is used for the numerical
  integration.
  %So shares the limitations of this numerical method. %\code{dkde} is a
  %shortcut for \code{kde(x=fhat$x, h=fhat$h,
  %eval.points=fhat$x)$estimate}.
  \code{rkde} uses
  Silverman (1986)'s method to generate a random sample from a KDE.
}

\references{ Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman \& Hall/CRC. London.}

\examples{
set.seed(8192)
x <- rnorm.mixt(n=10000, mus=0, sigmas=1, props=1)
fhat <- kde(x=x, h=hpi(x))
p1 <- pkde(fhat=fhat, q=c(-1, 0, 0.5))
qkde(fhat=fhat, p=p1)     ## should be close to c(-1, 0, 0.5)

x1 <- rkde(fhat, n=100)
plot(fhat)
fhat1 <- kde(x=x1, h=hpi(x1))
plot(fhat1, add=TRUE, col=2)
fhat2 <- dkde(x=x1, fhat=fhat1)
points(x1, fhat2, col=3)
}

\keyword{ smooth}

