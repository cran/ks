\name{Hamise.mixt, Hmise.mixt, Hamise.mixt.diag, Hmise.mixt.diag, amise.mixt, ise.mixt, mise.mixt}
\alias{Hmise.mixt}
\alias{Hamise.mixt}
\alias{Hmise.mixt.diag}
\alias{Hamise.mixt.diag}
\alias{hmise.mixt}
\alias{hamise.mixt}
\alias{ise.mixt}
\alias{amise.mixt}
\alias{mise.mixt}

\title{Squared error bandwidth matrix selectors for normal
  mixture densities}

\description{
  The global errors
  ISE (Integrated Squared Error), MISE (Mean Integrated Squared Error) and the 
  AMISE (Asymptotic Mean Integrated Squared Error) for 1- to 6-dimensional
  data.

  Normal mixture densities have closed form expressions for the MISE and
  AMISE. So in these cases, we can numerically minimise these criteria
  to find MISE- and AMISE-optimal matrices.
}


\usage{

Hamise.mixt(mus, Sigmas, props, samp, Hstart, deriv.order=0)
Hmise.mixt(mus, Sigmas, props, samp, Hstart, deriv.order=0)
Hamise.mixt.diag(mus, Sigmas, props, samp, Hstart, deriv.order=0)
Hmise.mixt.diag(mus, Sigmas, props, samp, Hstart, deriv.order=0)
hamise.mixt(mus, sigmas, props, samp, hstart, deriv.order=0)
hmise.mixt(mus, sigmas, props, samp, hstart, deriv.order=0)

ise.mixt(x, H, mus, Sigmas, props, h, sigmas, deriv.order=0)  
mise.mixt(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
amise.mixt(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors/vector of means}
  \item{sigmas, Sigmas}{vector of standard deviations/(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  \item{hstart, Hstart}{initial bandwidth (matrix), used in numerical
    optimisation}
  \item{deriv.order}{derivative order}
  \item{x}{matrix of data values}
  \item{h, H}{bandwidth (matrix)}
}


\value{
  -- Full MISE- or AMISE-optimal bandwidth matrix. Diagonal forms of
  these matrices are not available.

  -- ISE, MISE or AMISE value. 
}

\details{ 
  For normal mixture densities, ISE, MISE and AMISE 
  have exact formulas for all dimensions. See Chac\'on, Duong \& Wand (2008).

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k =
  \eqn{\left[\frac{4}{n(d+2r+2)}\right]^{2/(d+2r+4)}}{4/(n*(d + 2r +
	2))^(2/(d+ 2r+ 4))}, n = sample size, d = dimension of data, r=
  derivative order. The default for \code{hstart} is the square root of
  this expression.
}


\note{ISE is a random variable that depends on the data
  \code{x}. MISE and AMISE are non-random and don't
  depend on the data.}


\references{Chac\'on J.E., Duong, T. \& Wand, M.P. (2009). Asymptotics for
	general multivariate kernel density derivative
	estimators. \emph{Statistica Sinica}. Accepted. 
}

\examples{
## 1-d
mus <- c(0, 2)
sigmas <- c(1, sqrt(0.7))
props <- c(1/2, 1/2)
samp <- 1000
h <- hmise.mixt(mus, sigmas, props, samp, deriv.order=0)
x <- rnorm.mixt(n=samp, mus=mus, sigmas=sigmas, props=props)
ise.mixt(x=x, h=h, mus=mus, sigmas=sigmas, props=props)
mise.mixt(h=h, mus=mus, sigmas=sigmas, props=props, samp=samp)

## 2-d 
mus <- rbind(c(0,0), c(2,2))
Sigma <- matrix(c(1, 0.7, 0.7, 1), nr=2, nc=2) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
samp <- 100
H <- Hamise.mixt(mus, Sigmas, props, samp, deriv.order=2)
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)
ise.mixt(x=x, H=H, mus=mus, Sigmas=Sigmas, props=props, deriv.order=2)
amise.mixt(H=H, mus=mus, Sigmas=Sigmas, props=props, samp=samp, deriv.order=2)

}

\keyword{smooth}

