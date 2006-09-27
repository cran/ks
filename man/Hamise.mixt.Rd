\name{Hamise.mixt, Hmise.mixt}
\alias{Hmise.mixt}
\alias{Hamise.mixt}


\title{MISE- and AMISE-optimal bandwidth matrix selectors for normal
  mixture densities}
\description{
  Normal mixture densities have closed form expressions for the MISE and
  AMISE. So in these cases, we can numerically minimise these criteria
 to find MISE- and AMISE-optimal matrices. }
\usage{
Hmise.mixt(mus, Sigmas, props, samp, Hstart)
Hamise.mixt(mus, Sigmas, props, samp, Hstart)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
}


\value{Full MISE- or AMISE-optimal bandwidth matrix. Please note that
  diagonal forms of these matrices are not available.}

\details{ For normal mixture densities, the MISE and AMISE have exact
  formulas. See Wand \& Jones (1995).

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data.
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\examples{
mus <- rbind(c(0,0,0), c(2,2,2))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
samp <- 1000
Hmise.mixt(mus, Sigmas, props, samp)
Hamise.mixt(mus, Sigmas, props, samp)


}

\keyword{ smooth}

