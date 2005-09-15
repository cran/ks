\name{compare, compare.kda.diag.cv, compare.kda.cv, compare.pda.cv}

\alias{compare}
\alias{compare.kda.diag.cv}
\alias{compare.kda.cv}
\alias{compare.pda.cv}

\title{Comparisons for kernel and parametric discriminant analysis}
\description{
  Comparisons for kernel and parametric discriminant analysis.
}
\usage{
compare(x.group, est.group)
compare.kda.cv(x, x.group, bw="plugin", prior.prob=NULL, Hstart, ...)
compare.kda.diag.cv(x, x.group, bw="plugin", prior.prob=NULL, ...)
compare.pda.cv(x, x.group, type="quad", prior.prob=NULL)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{est.group}{vector of estimated group labels}
  \item{bw}{\code{"plugin"} = plug-in, \code{"lscv"} = LSCV, \code{"scv"} = SCV}
  \item{Hstart}{(stacked) matrix of intitial bandwidth matrices}
  \item{prior.prob}{vector of prior probabilities}
  \item{type}{\code{"line"} = linear discriminant, \code{"quad"} =
    quadratic discriminant}
  \item{...}{other optional parameters for bandwidth selection, see
    \code{\link{Hpi}}, \code{\link{Hlscv}}, \code{\link{Hscv}}} 
}

\value{ 
  The function \code{compare} creates a comparison between the true
  group labels \code{x.group} and the estimated ones \code{est.group}.
  It returns a list with fields
  \item{cross}{cross-classification table with the rows
    indicating the true group and the columns the estimated group}
  \item{error}{misclassification rate (MR)}
    
  In the case where we have test data that is independent of the
  training data, then  
    \deqn{\textrm{MR} = \frac{\textrm{number of points wrongly
	  classified}}{\textrm{total number of
	  points}}}{MR = (number of points wrongly classified) / (total number of
      points)}
    
  In the case where we don't have independent test data e.g.
  we are classifying the
  training data set, then the cross validated estimate is more
  appropriate.  See Silverman (1986). These are implemented as for
  kernel discriminant analysis as \code{compare.kda.cv} (full bandwidth
  selectors) and \code{compare.kda.diag.cv} (for diagonal bandwidth
  selectors), and \code{compare.pda.cv} for parametric discriminant analysis.  
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
  Otherwise the default is to use the sample proportions as 
  estimates of the prior probabilities.

  The parametric discriminant analysers use the code from the
  \code{MASS} library namely \code{lda} and \code{qda} for linear and
  quadratic discriminants.
}

\seealso{  \code{\link{kda.kde}}, \code{\link{pda.pde}}
}

\examples{
### bivariate example - restricted iris dataset  
library(MASS)
data(iris)
ir <- iris[,c(1,2)]
ir.gr <- iris[,5]

compare.kda.cv(ir, ir.gr, bw="plug-in", pilot="samse")
compare.pda.cv(ir, ir.gr, type="line")
compare.pda.cv(ir, ir.gr, type="quad")
}

\keyword{ smooth }
