\name{kda}
\alias{Hkda}
\alias{Hkda.diag}
\alias{kda}
\alias{hkda}
\alias{compare}
\alias{compare.kda.diag.cv}
\alias{compare.kda.cv}

\title{Kernel discriminant analysis for multivariate data}
\description{
  Kernel discriminant analysis for 1- to 6-dimensional data.
}
\usage{
Hkda(x, x.group, Hstart, bw="plugin", nstage=2, pilot="samse",
     pre="sphere", binned=FALSE, bgridsize)
Hkda.diag(x, x.group, bw="plugin", nstage=2, pilot="samse", 
     pre="sphere", binned=FALSE, bgridsize)
hkda(x, x.group, bw="plugin", nstage=2, binned=TRUE, bgridsize)

kda(x, x.group, Hs, hs, y, prior.prob=NULL)
compare(x.group, est.group, by.group=FALSE)
compare.kda.cv(x, x.group, bw="plugin", prior.prob=NULL, Hstart,
    by.group=FALSE, trace=FALSE, binned=FALSE, bgridsize,
    recompute=FALSE, ...)
compare.kda.diag.cv(x, x.group, bw="plugin", prior.prob=NULL,
    by.group=FALSE, trace=FALSE, binned=FALSE, bgridsize,
    recompute=FALSE, ...)
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
  \item{pilot}{pilot selector: see \code{\link{Hpi}}, \code{\link{Hscv}}}
  \item{pre}{\code{"scale"} = \code{\link{pre.scale}}, \code{"sphere"} = \code{\link{pre.sphere}}}
  \item{Hstart}{(stacked) matrix of initial bandwidth matrices, used in
    numerical optimisation}
  \item{binned}{flag for binned kernel estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{est.group}{vector of estimated group labels}
  \item{by.group}{flag to give results also within each group}
  \item{trace}{flag for printing messages in command line to
    trace the execution}
  \item{recompute}{flag for recomputing the bandwidth matrix after
    excluding the i-th data item}
  \item{...}{other optional parameters for bandwidth selection, see
    \code{\link{Hpi}}, \code{\link{Hlscv}}, \code{\link{Hscv}}} 
}

\value{
  --The result from \code{Hkda} and \code{Hkda.diag} is a stacked matrix
  of bandwidth matrices, one for each training data group. The result
  from \code{hkda} is a vector of bandwidths, one for each training data
  group.  

  --The result from \code{kda} is a vector of group labels
  estimated via the kernel discriminant rule.  If the test data \code{y} are
  given then these are classified. Otherwise the training data \code{x}
  are classified.

  --The \code{compare} functions create a comparison between the true
  group labels \code{x.group} and the estimated ones. 
  It returns a list with fields
  \item{cross}{cross-classification table with the rows
    indicating the true group and the columns the estimated group}
  \item{error}{misclassification rate (MR)}
    
  In the case where we have test data that is independent of the
  training data, \code{compare} computes MR = (number of points wrongly classified)/
  (total number of points).
    
  In the case where we don't have independent test data e.g.
  we are classifying the training data set itself, then the cross
  validated estimate of MR is more appropriate. These
  are implemented as \code{compare.kda.cv} (full bandwidth
  selectors) and \code{compare.kda.diag.cv} (for diagonal bandwidth
  selectors). These functions are only available for d > 1.

  If \code{by.group=FALSE} then only the total MR rate is given. If it
  is set to TRUE, then the MR rates for each class are also given
  (estimated number in group divided by true number).
}

\references{
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York
}

\details{
  --The values that valid for \code{bw} are \code{"plugin", "lscv"} and
  \code{"scv"} for
  \code{Hkda}.  These in turn call \code{\link{Hpi}},
  \code{\link{Hlscv}} and \code{\link{Hscv}}. For plugin selectors, all
  of \code{nstage}, \code{pilot} and \code{pre} need to be set. For SCV
  selectors, currently \code{nstage=1} always but \code{pilot} and \code{pre}
  need to be set.  For LSCV selectors, none of them are required. 
  \code{Hkda.diag} makes analogous calls to diagonal selectors.

  --If you have prior probabilities then set \code{prior.prob} to these.
  Otherwise \code{prior.prob=NULL} is the default i.e. use the sample
  proportions as estimates of the prior probabilities.
  
  If \code{trace=TRUE}, a message is printed in the command line
  indicating that it's processing the i-th data item:
  cross-validated estimates may take a long time to execute.
}

\examples{ 
### See examples in ? plot.kda.kde

### univariate example -- independent test data
x <- c(rnorm.mixt(n=100, mus=1, sigmas=1, props=1),
       rnorm.mixt(n=100, mus=-1, sigmas=1, props=1))
x.gr <- rep(c(1,2), times=c(100,100))
y <- c(rnorm.mixt(n=100, mus=1, sigmas=1, props=1),
       rnorm.mixt(n=100, mus=-1, sigmas=1, props=1))
kda.gr <- kda(x, x.gr, hs=sqrt(c(0.09, 0.09)), y=y)
compare(x.gr, kda.gr)
compare(x.gr, kda.gr, by.group=TRUE) 

### bivariate example - restricted iris dataset, dependent test data
library(MASS)
data(iris)
ir <- iris[,c(1,2)]
ir.gr <- iris[,5]
compare.kda.cv(ir, ir.gr, bw="plug-in", pilot="samse") } 

\keyword{ smooth }
