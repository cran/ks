\name{Hmise.mixt, Hamise.mixt}
\alias{Hmise.mixt}
\alias{Hamise.mixt}


\title{MISE- and AMISE-optimal bandwidth matrix selectors for normal
  mixture densities}
\description{
  For normal mixture densities, we have a closed form for the MISE and
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
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
samp <- 100
Hmise.mixt(mus, Sigmas, props, samp)
Hamise.mixt(mus, Sigmas, props, samp)
}

\keyword{ smooth}

