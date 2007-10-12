\name{compare, compare.kda.diag.cv, compare.kda.cv}

\alias{compare}
\alias{compare.kda.diag.cv}
\alias{compare.kda.cv}

\title{Comparisons for kernel discriminant analysis}
\description{
  Comparisons for kernel discriminant analysis.
}
\usage{
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
  \item{est.group}{vector of estimated group labels}
  \item{bw}{\code{"plugin"} = plug-in, \code{"lscv"} = LSCV, \code{"scv"} = SCV}
  \item{Hstart}{(stacked) matrix of initial bandwidth matrices}
  \item{prior.prob}{vector of prior probabilities}
  \item{by.group}{if set to TRUE then results are given for each group
    as well}
  \item{trace}{if set to TRUE then message is printed in command line to
    trace the execution}
  \item{binned}{flag for binned kernel estimation}
  \item{bgridsize}{vector of binning grid sizes - only required if
    \code{binned=TRUE}}
  \item{recompute}{flag for recomputing the bandwidth matrix after
    excluding the i-th data item}
  \item{...}{other optional parameters for bandwidth selection, see
    \code{\link{Hpi}}, \code{\link{Hlscv}}, \code{\link{Hscv}}} 
}

\value{ 
  The functions create a comparison between the true
  group labels \code{x.group} and the estimated ones. 
  It returns a list with fields
  \item{cross}{cross-classification table with the rows
    indicating the true group and the columns the estimated group}
  \item{error}{misclassification rate (MR)}
    
  In the case where we have test data that is independent of the
  training data, \code{compare} computes    
    \deqn{\textrm{MR} = \frac{\textrm{number of points wrongly
	  classified}}{\textrm{total number of
	  points}}.}{MR = (number of points wrongly classified) / (total number of
      points).}
    
  In the case where we don't have independent test data e.g.
  we are classifying the training data set itself, then the cross
  validated estimate of MR is more appropriate.  See Silverman (1986). These
  are implemented as \code{compare.kda.cv} (full bandwidth
  selectors) and \code{compare.kda.diag.cv} (for diagonal bandwidth
  selectors). These functions are only available for d > 1.

  If \code{by.group=FALSE} then only the total MR rate is given. If it
  is set to TRUE, then the MR rates for each class are also given
  (estimated number in group divided by true number).

}

\references{

  Silverman, B. W. (1986) \emph{Data Analysis for Statistics and Data
    Analysis}. Chapman \& Hall. London.
  
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York

  Venables, W.N. & Ripley, B.D. (1997) \emph{Modern Applied Statistics with
    S-PLUS}. Springer-Verlag. New York.   
}

\details{
  If you have prior probabilities then set \code{prior.prob} to these.
  Otherwise \code{prior.prob=NULL} is the default i.e. use
  the sample proportions as 
  estimates of the prior probabilities.

  If \code{trace=TRUE}, a message is printed in the command line
  indicating that it's processing the i-th data item. This can be
  helpful since the cross-validated estimates may take a long time to
  execute completely.
}

\seealso{  \code{\link{kda.kde}}}
}

\examples{
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

compare.kda.cv(ir, ir.gr, bw="plug-in", pilot="samse")
compare.kda.diag.cv(ir, ir.gr, bw="plug-in", pilot="samse",
   pre="scaled", binned=TRUE)

\dontrun{
compare.kda.cv(ir, ir.gr, bw="plug-in", pilot="samse", recompute=TRUE)
}

}

\keyword{ smooth }
