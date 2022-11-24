\name{Hnm}
\alias{Hnm}
\alias{hnm}
\alias{Hnm.diag}
\title{Normal mixture bandwidth}

\description{
Normal mixture bandwidth.  
}

\usage{
Hnm(x, deriv.order=0, G=1:9, subset.ind, mise.flag=FALSE, verbose, ...)
Hnm.diag(x, deriv.order=0, G=1:9, subset.ind, mise.flag=FALSE, verbose, ...)
hnm(x, deriv.order=0, G=1:9, subset.ind, mise.flag=FALSE, verbose, ... )
}

\arguments{
  \item{x}{vector/matrix of data values}
  \item{deriv.order}{derivative order}
  \item{G}{range of number of mixture components}
  \item{subset.ind}{index vector of subset of \code{x} for fitting}
  \item{mise.flag}{flag to use MISE or AMISE minimisation. Default is FALSE.}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{...}{other parameters for \code{Mclust}}
}

\value{
  Normal mixture bandwidth. If \code{mise=TRUE} then the minimal MISE value is returned too. 
}

\details{
  The normal mixture fit is provided by the \code{Mclust} function in
  the \pkg{mclust} package. \code{Hnm} is then \code{Hmise.mixt} (if
  \code{mise.flag=TRUE}) or \code{Hamise.mixt} (if \code{mise.flag=FALSE}) with
  these fitted normal mixture parameters. Likewise for 
  \code{Hnm.diag}, \code{hnm}.   
}

\seealso{ \link{Hmise.mixt}, \link{Hamise.mixt} }

\references{Cwik, J. & Koronacki, J.  (1997). A combined
  adaptive-mixtures/plug-in estimator of multivariate probability
  densities. \emph{Computational Statistics and Data Analysis},
  \bold{26}, 199-218.   
}

\examples{
data(unicef)
Hnm(unicef) 
}

\keyword{smooth}

