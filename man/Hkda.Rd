\name{Hkda, Hkda.diag}
\alias{Hkda}
\alias{Hkda.diag}

\title{Bandwidth matrix selectors for kernel discriminant analysis for
  bivariate data}
\description{
  Bandwidth matrices for kernel discriminant analysis for bivariate data.
}
\usage{
Hkda(x, x.group, Hstart, bw="plugin", nstage=2, pilot="samse",pre="sphere")
Hkda.diag(x, x.group, bw="plugin", nstage=2, pilot="samse", pre="sphere")
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{bw}{bandwidth: \code{"plugin"} = plug-in, \code{"lscv"} = LSCV, 
  \code{"scv"} = SCV} 
 \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"amse"} = AMSE-optimal pilot bandwidths,
    \code{"samse"} = single SAMSE-optimal pilot bandwidth}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} =
    pre-sphering}
  \item{Hstart}{(stacked) matrix of initial bandwidth matrices, used in
    numerical optimisation}
}

\value{
  Stacked matrix of bandwidth matrices for each training data group.
}
\references{
  Duong, T. (2004) \emph{Bandwidth Matrices for Multivariate Kernel Density 
  Estimation.} Ph.D. Thesis. University of Western Australia. 

  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.
}

\details{
  The values that valid for \code{bw} are \code{"plugin", "lscv"} and
  \code{"scv"} for
  \code{Hkda}.  These in turn call \code{\link{Hpi}},
  \code{\link{Hlscv}} and \code{\link{Hscv}}. For plugin selectors, all
  of \code{nstage}, \code{pilot} and \code{pre} need to be set. For SCV
  selectors, currently \code{nstage} is always programmed to be one but
  the other two need to be set.  For LSCV selectors, none of them are required. 


  For \code{Hkda.diag}, only \code{"plugin"} or \code{"lscv"} are valid
  which in turn call  \code{\link{Hpi.diag}} 
  and \code{\link{Hlscv.diag}}. Again, \code{nstage}, \code{pilot} and
  \code{pre} are available for \code{Hpi.diag} but not required for
  \code{Hlscv.diag}.
}

\seealso{ \code{\link{kda.kde}}, \code{\link{Hpi}}, \code{\link{Hpi.diag}},
    \code{\link{Hlscv}}, \code{\link{Hlscv.diag}}, \code{\link{Hscv}}
}

\examples{
library(MASS)
data(iris)
iris.mat <- rbind(iris[,,1], iris[,,2], iris[,,3])
ir <- iris.mat[,c(1,2)]
ir.gr <- iris.mat[,5]

Hkda(ir, ir.gr, bw="scv", pre="scale")
Hkda.diag(ir, ir.gr, bw="plugin", pre="scale")
}

\keyword{ smooth }
