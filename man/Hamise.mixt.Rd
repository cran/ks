\name{Hamise.mixt, Hmise.mixt, amise.mixt, ise.mixt, mise.mixt}
\alias{Hmise.mixt}
\alias{Hamise.mixt}
\alias{ise.mixt}
\alias{mise.mixt}
\alias{amise.mixt}

\title{MISE- and AMISE-optimal bandwidth matrix selectors for normal
  mixture densities}
\description{
  Normal mixture densities have closed form expressions for the MISE and
  AMISE. So in these cases, we can numerically minimise these criteria
  to find MISE- and AMISE-optimal matrices.
  
  The global errors
  ISE (Integrated Squared Error), MISE (Mean Integrated Squared Error)
  of kernel density estimates
  for normal densities, for 2- to 6-dimensional data, and
  AMISE (Asymptotic Mean Integrated Squared Error) for 2-dimensional
  data.}


\usage{
Hmise.mixt(mus, Sigmas, props, samp, Hstart)
Hamise.mixt(mus, Sigmas, props, samp, Hstart)

ise.mixt(x, H, mus, Sigmas, props)  
mise.mixt(H, mus, Sigmas, props, samp)
amise.mixt(H, mus, Sigmas, props, samp)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
}


\value{
  -- Full MISE- or AMISE-optimal bandwidth matrix. Diagonal forms of
  these matrices are not available.

  -- ISE, MISE or AMISE value.
}

\details{ 
  For normal mixture densities, ISE and MISE 
  have exact formulas for all dimensions, and AMISE has an exact form for
  2 dimensions. See Wand \& Jones (1995).

  If \code{Hstart} is not given then it defaults to
  \code{k*var(x)} where k = \eqn{\left[\frac{4}{n(d+2)}\right]^{2/(d+4)}}{4/(n*(d + 2))^(2/(d+ 4))}, n = sample size, d = dimension of data.
}


\note{ISE is a random variable that depends on the data
  \code{x}. MISE and AMISE are non-random and don't
  depend on the data.}


\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\examples{
mus <- rbind(c(0,0,0), c(2,2,2))
Sigma <- matrix(c(1, 0.7, 0.7, 0.7, 1, 0.7, 0.7, 0.7, 1), nr=3, nc=3) 
Sigmas <- rbind(Sigma, Sigma)
props <- c(1/2, 1/2)
samp <- 1000
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)

H1 <- Hmise.mixt(mus, Sigmas, props, samp)
H2 <- Hamise.mixt(mus, Sigmas, props, samp)

ise.mixt(x, H2, mus, Sigmas, props)
mise.mixt(H2, mus, Sigmas, props, samp)
}

\keyword{ smooth}

