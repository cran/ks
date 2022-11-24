\name{rkde}
\alias{dkde}
\alias{pkde}
\alias{qkde}
\alias{rkde}
\title{Derived quantities from kernel density estimates}

\description{
   Derived quantities from kernel density estimates.
}

\usage{
 dkde(x, fhat)
 pkde(q, fhat)
 qkde(p, fhat)
 rkde(n, fhat, positive=FALSE) 
}

\arguments{
   \item{x,q}{vector of quantiles}
   \item{p}{vector of probabilities}
   \item{n}{number of observations}
   \item{positive}{flag to compute KDE on the positive real
   line. Default is FALSE.}  
   \item{fhat}{kernel density estimate, object of class \code{kde}}
}

\value{
   For the 1-d kernel density estimate \code{fhat},
   \code{pkde} computes the cumulative probability for the quantile
   \code{q}, \code{qkde} computes the quantile corresponding to the probability
   \code{p}.

   For any kernel density estimate, \code{dkde} computes the density value at
   \code{x} (it is an alias for \code{predict.kde}), \code{rkde}
   computes a random sample of size \code{n}.   
}

\details{ 
  \code{pkde} uses the trapezoidal rule for the numerical
  integration. \code{rkde} uses
  Silverman (1986)'s method to generate a random sample from a KDE. 
}

\references{ Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman & Hall/CRC. London.}

\examples{
set.seed(8192)
x <- rnorm.mixt(n=10000, mus=0, sigmas=1, props=1)
fhat <- kde(x=x)
p1 <- pkde(fhat=fhat, q=c(-1, 0, 0.5))
qkde(fhat=fhat, p=p1)    
y <- rkde(fhat=fhat, n=100)

x <- rmvnorm.mixt(n=10000, mus=c(0,0), Sigmas=invvech(c(1,0.8,1)))
fhat <- kde(x=x)
y <- rkde(fhat=fhat, n=1000)
fhaty <- kde(x=y)
plot(fhat, col=1)
plot(fhaty, add=TRUE, col=2)
}

\keyword{smooth}

