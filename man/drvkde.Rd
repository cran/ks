\name{drvkde}
\alias{drvkde}
\title{Kernel density derivative estimation}

\description{
  Compute kernel density derivative
  estimates and standard errors for multivariate data. 
}

\usage{
drvkde(x, drv, bandwidth, gridsize, range.x, binned=FALSE, se=TRUE)
}

\arguments{
  \item{x}{data matrix or matrix of binning counts} 
  \item{drv}{vector of derivative indices}
  \item{bandwidth}{vector of bandwidths}
  \item{gridsize}{vector of grid sizes}
  \item{range.x}{list of vector of ranges for x}
  \item{binned}{TRUE if x is binned counts or FALSE if x is
    data matrix} 
  \item{se}{flag for computing the standard error of kernel
    estimate}
 %%\item{estimate.positive}{flag for ensuring all estimate values are positive}
}

\value{
  Returns a list with fields
  
  \code{x.grid}  - grid points \cr
  \code{est} - kernel estimate of partial derivative of density function
  indicated by \code{drv}\cr
  \code{se} -  estimate of standard error of \code{est} (if \code{se=TRUE}).
}

\details{The estimates and standard errors are computed over a grid of binned counts
  \code{x.grid}. If the binned counts are not supplied then they are computed
  inside this function.

  If \code{gridsize} and \code{range.x} are not supplied, they are
  computed inside this function.
}


\examples{
## univariate
x <- rnorm(100)
fhat <- drvkde(x=x, drv=0, bandwidth=0.1)    ## KDE of f
fhat1 <- drvkde(x=x, drv=1, bandwidth=0.1)   ## KDE of df/dx
}

\references{
     Wand, M.P. and Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman
	 \& Hall/CRC, London.
}
   
\author{M.P. Wand}


\keyword{smooth} 

