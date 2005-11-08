\name{kda.kde, pda.pde}
\alias{kda.kde}
\alias{pda.pde}

\title{Density estimates for kernel and parametric discriminant analysis
   for multivariate data}
\description{
  Density estimates for kernel and parametric discriminant analysis 
   for 2- to 6-dimensional data.
}
\usage{
kda.kde(x, x.group, Hs, gridsize, supp=3.7, eval.points=NULL)
pda.pde(x, x.group, gridsize, type="quad", xlim, ylim, zlim)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{Hs}{(stacked) matrix of bandwidth matrices}
  \item{gridsize}{vector of number of grid points}
  \item{supp}{effective support for standard normal is [\code{-supp, supp}]}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{type}{\code{"line"} = linear discriminant, \code{"quad"} =
    quadratic discriminant}
  \item{xlim, ylim, zlim}{x-axis, y-axis, z-axis limits}
}
  

\value{
  Density estimate for discriminant analysis is an object of
  class \code{dade} which is a
  list with 6 fields
  \item{x}{data points - same as input}
  \item{eval.points}{points that density estimate is evaluated at}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{H}{bandwidth matrices}
  \item{prior.prob}{sample proportions of each group}
  \item{type}{one of \code{"kernel", "linear", "quadratic"} indicating
      the type of discriminant analyser used.}
}

\details{
  The kernel density estimate is based on \code{\link{kde}}.
  
  If \code{gridsize} is not set to a specific value, then it
  defaults to 100 grid points in each co-ordinate direction
  i.e. \code{rep(100, d)}. Not required
  to be set if specifying \code{eval.points}.

  If \code{eval.points} is not specified, then the
  density estimate is automatically computed over a grid whose
  resolution is controlled by \code{gridsize} (a grid is
  required for plotting).

  The parametric discriminant analysers use the code from the
  \code{MASS} library namely \code{lda} and \code{qda} for linear and
  quadratic discriminants. 

  If \code{xlim} and \code{ylim} are not specified then they default to
  be 10 \% bigger than the range of the data values.
}

\references{
  Mardia, K.V., Kent, J.T. \& Bibby J.M. (1979) \emph{Multivariate
    Analysis}. Academic Press. London.
  
  Simonoff, J. S., (1996) \emph{Smoothing Methods in Statistics},
  Springer-Verlag. New York.
  
  Venables, W.N. & Ripley, B.D. (1997) \emph{Modern Applied Statistics with
    S-PLUS}, Springer-Verlag. New York. }
}

\seealso{\code{\link{plot.dade}}, \code{\link{pda}}, \code{\link{kda}},
  \code{\link{kde}}
}

\examples{
library(MASS)
data(iris)
ir <- iris[,c(1,2)]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, H, ir)
fhat <- kda.kde(ir, ir.gr, H)

qda.gr <- pda(ir, ir.gr, ir, type="quad")
qda.fhat <- pda.pde(ir, ir.gr)
}

\keyword{smooth}
