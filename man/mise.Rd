\name{ise, mise, amise}
\alias{ise.mixt}
\alias{mise.mixt}
\alias{amise.mixt}


\title{ISE, MISE and AMISE of kernel density estimates for normal mixture densities}
\description{
  The global errors
  ISE (Integrated Squared Error), MISE (Mean Integrated Squared Error)
  of kernel density estimates
  for normal densities, for 2- to 6-dimensional data, and
and AMISE (Asymptotic Mean Integrated Squared Error) for 2-dimensional data.}

\usage{
ise.mixt(x, H, mus, Sigmas, props)  
mise.mixt(H, mus, Sigmas, props, samp)
amise.mixt(H, mus, Sigmas, props, samp)
}

\arguments{
  \item{x}{matrix of data values}
  \item{H}{bandwidth matrix}
  \item{mus}{(stacked) matrix of mean vectors}
  \item{Sigmas}{(stacked) matrix of variance matrices}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  }
  

\value{ISE, MISE or AMISE value.}

\details{ For normal mixture densities, ISE and MISE 
  have exact formulas for all dimensions, and AMISE has an exact form for
  2 dimensions. See Wand \& Jones (1995).  
}

\references{ Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}.
  Chapman \& Hall. London. } 

\note{Remember that ISE is a random variable that depends on the data
  \code{x}; and that MISE and AMISE are non-random and don't
  depend on the data.}

\examples{
samp <- 50
d <- 3
mus <- rbind(rep(0,d), rep(1,d))
Sigmas <- 0.25*rbind(diag(d), diag(d))
props <- c(2/3, 1/3)
x <- rmvnorm.mixt(samp, mus, Sigmas, props)
H <- Hpi(x)
ise.mixt(x, H, mus, Sigmas, props)
mise.mixt(H, mus, Sigmas, props, samp)
}

\keyword{ smooth}
