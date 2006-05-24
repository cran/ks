\name{kda, pda, Hkda, Hkda.diag, kda.kde, pda.pde}
\alias{Hkda}
\alias{Hkda.diag}
\alias{kda}
\alias{pda}
\alias{kda.kde}
\alias{pda.pde}

\title{Kernel and parametric discriminant analysis}
\description{
  Kernel and parametric discriminant analysis.
}
\usage{
Hkda(x, x.group, Hstart, bw="plugin", nstage=2, pilot="samse",
     pre="sphere", binned=FALSE)
Hkda.diag(x, x.group, bw="plugin", nstage=2, pilot="samse", 
     pre="sphere", binned=FALSE)

kda(x, x.group, Hs, y, prior.prob=NULL)
pda(x, x.group, y, prior.prob=NULL, type="quad")

kda.kde(x, x.group, Hs, gridsize, supp=3.7, eval.points=NULL)
pda.pde(x, x.group, gridsize, type="quad", xlim, ylim, zlim)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{y}{matrix of test data} 
  \item{Hs}{(stacked) matrix of bandwidth matrices}
  \item{prior.prob}{vector of prior probabilities}
  \item{type}{\code{"line"} = linear discriminant, \code{"quad"} =
    quadratic discriminant}

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

  \item{gridsize}{vector of number of grid points}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{xlim, ylim, zlim}{x-axis, y-axis, z-axis limits (used only for plotting)}
}

\value{
  -- The result from \code{Hkda} and \code{Hkda.diag} is a stacked matrix
  of bandwidth matrices for each training data group. This is then
  suitable to passed as the \code{Hs} argument in \code{kda}.

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
  
  -- The result from \code{kda} and \code{pda} is a vector of group labels
  estimated via a discriminant (or classification) 
  rule.  If the test data \code{y} are given then these are
  classified. Otherwise the training data \code{x} are classified.

  -- The result from \code{kda.kde} and \code{pda.pde} is a density estimate
  for discriminant analysis is an object of  class \code{dade} which is a
  list with 6 fields
  \item{x}{data points - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrices}
  \item{prior.prob}{sample proportions of each group}
  \item{type}{one of \code{"kernel", "linear", "quadratic"} indicating
      the type of discriminant analyser used.}
  
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

  The linear and quadratic discriminant analysers are based on
  \code{lda} and \code{qda} from the \code{MASS} library.

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

  -- The kernel density estimate is based on \code{\link{kde}}.
  
  If \code{eval.points=NULL} (default) then the
  density estimate is automatically computed over a grid whose
  resolution is controlled by \code{gridsize} (default is
  100 in each co-ordinate direction).
  
  If \code{xlim} and \code{ylim} are not specified then they default to
  be 10\% bigger than the range of the data values.
  
}

\seealso{
  \code{\link{compare}}, 
  \code{\link{compare.kda.cv}},
  \code{\link{compare.pda.cv}} 
}

\examples{

### bivariate example - restricted iris dataset  
library(MASS)
data(iris)
ir <- iris[,1:2]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, H, ir)
lda.gr <- pda(ir, ir.gr, ir, type="line")
qda.gr <- pda(ir, ir.gr, ir, type="quad")

\dontrun{
### multivariate example - full iris dataset
ir <- iris[,1:4]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, H, ir)
lda.gr <- pda(ir, ir.gr, ir, type="line")
qda.gr <- pda(ir, ir.gr, ir, type="quad")
}
}
\keyword{ smooth }
