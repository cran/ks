\name{Hamise.mixt}
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
  data. Normal mixture densities have closed form expressions for the MISE and
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

amise.mixt(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
ise.mixt(x, H, mus, Sigmas, props, h, sigmas, deriv.order=0, 
         binned=FALSE, bgridsize)  
mise.mixt(H, mus, Sigmas, props, samp, h, sigmas, deriv.order=0)
}

\arguments{
  \item{mus}{(stacked) matrix of mean vectors (>1-d), vector of means (1-d)}
  \item{Sigmas,sigmas}{(stacked) matrix of variance matrices (>1-d), vector of standard deviations (1-d)}
  \item{props}{vector of mixing proportions}
  \item{samp}{sample size}
  \item{Hstart,hstart}{initial bandwidth (matrix), used in numerical optimisation}
  \item{deriv.order}{derivative order}
  \item{x}{matrix of data values}
  \item{H,h}{bandwidth (matrix)}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
}


\value{
  Full MISE- or AMISE-optimal bandwidth matrix. ISE, MISE or AMISE value. 
}

\details{
  ISE is a random variable that depends on the data
  \code{x}. MISE and AMISE are non-random and don't
  depend on the data. For normal mixture densities, ISE, MISE and AMISE 
  have exact formulas for all dimensions. See Chacon, Duong & Wand (2011).
}


\references{Chacon J.E., Duong, T. & Wand, M.P. (2011). Asymptotics for
	general multivariate kernel density derivative
	estimators. \emph{Statistica Sinica}. \bold{21}, 807-840.  
}

\examples{
mus <- rbind(c(0,0), c(2,2))
Sigmas <- rbind(invvech(c(1, 0.7, 1)), invvech(c(1, 0.7, 1)))
props <- c(1/2, 1/2)
samp <- 100
H <- Hamise.mixt(mus, Sigmas, props, samp, deriv.order=2)
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)
ise.mixt(x=x, H=H, mus=mus, Sigmas=Sigmas, props=props, deriv.order=2)
}

\keyword{smooth}

