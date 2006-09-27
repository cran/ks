\name{kda, Hkda, Hkda.diag}
\alias{Hkda}
\alias{Hkda.diag}
\alias{kda}

\title{Kernel discriminant analysis for multivariate data}
\description{
  Kernel discriminant analysis for 1- to 6-dimensional data.
}
\usage{
kda(x, x.group, Hs, hs, y, prior.prob=NULL)

Hkda(x, x.group, Hstart, bw="plugin", nstage=2, pilot="samse",
     pre="sphere", binned=FALSE)
Hkda.diag(x, x.group, bw="plugin", nstage=2, pilot="samse", 
     pre="sphere", binned=FALSE)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{y}{matrix of test data} 
  \item{Hs}{(stacked) matrix of bandwidth matrices}
  \item{hs}{vector of scalar bandwidths}
  \item{prior.prob}{vector of prior probabilities}
  \item{bw}{bandwidth: \code{"plugin"} = plug-in, \code{"lscv"} = LSCV, 
    \code{"scv"} = SCV} 
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"amse"} = AMSE pilot bandwidths,
    \code{"samse"} = single SAMSE pilot bandwidth}
  \item{pre}{\code{"scale"} = pre-scaling, \code{"sphere"} =
    pre-sphering}
  \item{Hstart}{(stacked) matrix of initial bandwidth matrices, used in
    numerical optimisation}
  \item{binned}{if TRUE (default) used binned estimation otherwise use
    exact kernel density estimation}
}

\value{
  -- The result from \code{kda} is a vector of group labels
  estimated via the kernel discriminant rule.  If the test data \code{y} are
  given then these are classified. Otherwise the training data \code{x}
  are classified.

  -- The result from \code{Hkda} and \code{Hkda.diag} is a stacked matrix
  of bandwidth matrices for each training data group. 

  The values that valid for \code{bw} are \code{"plugin", "lscv"} and
  \code{"scv"} for \code{Hkda}.  These in turn call \code{\link{Hpi}},
  \code{\link{Hlscv}} and \code{\link{Hscv}}. For plugin selectors, all
  of \code{nstage}, \code{pilot} and \code{pre} need to be set. For SCV
  selectors, currently \code{nstage=1} always but \code{pilot} and \code{pre}
  need to be set.  For LSCV selectors, none of them are required. 

  For \code{Hkda.diag}, options are \code{"plugin"} or \code{"lscv"}
  which in turn call respectively \code{\link{Hpi.diag}} 
  and \code{\link{Hlscv.diag}}. Again, \code{nstage}, \code{pilot} and
  \code{pre} are available for \code{Hpi.diag} but not required for
  \code{Hlscv.diag}.
  
  For details on the pre-transformations in \code{pre}, see
  \code{\link{pre.sphere}} and \code{\link{pre.scale}}.
  

}

\references{
  Mardia, K.V., Kent, J.T. \& Bibby J.M. (1979) \emph{Multivariate
    Analysis}. Academic Press. London.
 
  Silverman, B. W. (1986) \emph{Data Analysis for Statistics and Data
    Analysis}. Chapman \& Hall. London.
  
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York

  Venables, W.N. & Ripley, B.D. (1997) \emph{Modern Applied Statistics with
    S-PLUS}. Springer-Verlag. New York.   
}

\details{
  -- If you have prior probabilities then set \code{prior.prob} to these.
  Otherwise \code{prior.prob=NULL} is the default i.e. use the sample
  proportions as estimates of the prior probabilities.

  -- The values that valid for \code{bw} are \code{"plugin", "lscv"} and
  \code{"scv"} for
  \code{Hkda}.  These in turn call \code{\link{Hpi}},
  \code{\link{Hlscv}} and \code{\link{Hscv}}. For plugin selectors, all
  of \code{nstage}, \code{pilot} and \code{pre} need to be set. For SCV
  selectors, currently \code{nstage=1} always but \code{pilot} and \code{pre}
  need to be set.  For LSCV selectors, none of them are required. 

  For \code{Hkda.diag}, options are \code{"plugin"} or \code{"lscv"}
  which in turn call respectively \code{\link{Hpi.diag}} 
  and \code{\link{Hlscv.diag}}. Again, \code{nstage}, \code{pilot} and
  \code{pre} are available for \code{Hpi.diag} but not required for
  \code{Hlscv.diag}.
  
}

\seealso{
  \code{\link{compare}}, 
  \code{\link{compare.kda.cv}},
  \code{\link{kda.kde}}
}

\examples{
### univariate example -- independent test data
x <- c(rnorm.mixt(n=100, mus=1, sigmas=1, props=1),
       rnorm.mixt(n=100, mus=-1, sigmas=1, props=1))
x.gr <- rep(c(1,2), times=c(100,100))
y <- c(rnorm.mixt(n=100, mus=1, sigmas=1, props=1),
       rnorm.mixt(n=100, mus=-1, sigmas=1, props=1))

kda.gr <- kda(x, x.gr, hs=sqrt(c(0.09, 0.09)), y=y)


### bivariate example - restricted iris dataset  
library(MASS)
data(iris)
ir <- iris[,1:2]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, Hs=H, y=ir)


\dontrun{
### multivariate example - full iris dataset
ir <- iris[,1:4]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, Hs=H, y=ir)
}
}
\keyword{ smooth }
