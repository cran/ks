\name{ise, mise, amise}
\alias{ise.mixt}
\alias{iset.mixt}
\alias{mise.mixt}
\alias{amise.mixt}


\title{ISE, MISE and AMISE of kernel density estimates for normal and t
  mixture densities}
\description{
  The global errors
  ISE (Integrated Squared Error), MISE (Mean Integrated Squared Error)
  and AMISE (Asymptotic Mean Integrated Squared Error) of kernel density estimates
  for normal and t mixture densities.}

\usage{
ise.mixt(x, H, mus, Sigmas, props, lower, upper, gridsize=c(250,250),
         stepsize)
iset.mixt(x, H, mus, Sigmas, dfs, props, lower, upper, gridsize=c(250,250),
          stepsize)  
mise.mixt(H, mus, Sigmas, props, samp)
amise.mixt(H, mus, Sigmas, props, samp)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{dfs}{vector of degrees of freedom}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  \item{lower, upper}{vectors of lower, upper bounds for numerical integration}
  \item{gridsize}{vector of number of points in each dimension}
  \item{stepsize}{vector of step sizes in each dimension}
  }
  

\value{ISE, MISE or AMISE value.}

\details{ For normal mixture densities, the ISE, MISE and AMISE all
  have exact formulas.  See Wand \& Jones (1995).  For the t mixture
  densities, we resort to using numerical integration, using a simple
  Riemann sum.  A grid is set up and the function values are computed
  and then multiplied by the area of the grid element to give an
  approximation of the volume under the curve.  The resolution of the
  grid is given either by \code{gridsize} or \code{stepsize}. 

}
\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\note{Remember that ISE is a random variable that depends on the data
  \code{x}; and that MISE and AMISE are non-random and don't
  depend on the data.}

\examples{
samp <- 100
mus <- rbind(c(-3/2,0), c(3/2,0))
Sigmas <- rbind(diag(c(1/16, 1)), rbind(c(1/16, 1/18), c(1/18, 1/16)))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(samp, mus, Sigmas, props)
H <- Hpi(x)
ise.mixt(x, H, mus, Sigmas, props, stepsize=0.01)
mise.mixt(H, mus, Sigmas, props, samp)
amise.mixt(H, mus, Sigmas, props, samp)

dfs <- c(7,5)
x <- rmvt.mixt(samp, mus, Sigmas, dfs, props)
H <- Hpi(x)
iset.mixt(x, H, mus, Sigmas, dfs, props, lower=c(-5,-5), upper=c(5,5))
}

\keyword{ smooth}
