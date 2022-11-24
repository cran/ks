\name{histde}
\alias{histde}
\alias{predict.histde}
\title{Histogram density estimate}

\description{
  Histogram density estimate for 1- and 2-dimensional data.
}

\usage{
histde(x, binw, xmin, xmax, adj=0)

\method{predict}{histde}(object, ..., x)
}

\arguments{
  \item{x}{matrix of data values}
  \item{binw}{(vector) of binwidths}
  \item{xmin,xmax}{vector of minimum/maximum values for grid}
  \item{adj}{displacement of default anchor point, in percentage of 1
    bin}
  \item{object}{object of class \code{histde}}
  \item{...}{other parameters}
}
  
\value{
  A histogram density estimate is an object of class \code{histde} which is a
  list with fields: 
  \item{x}{data points - same as input}
  \item{eval.points}{vector or list of points at which the estimate is evaluated}
  \item{estimate}{density estimate at \code{eval.points}}
  \item{binw}{(vector of) bandwidths}
  \item{nbin}{(vector of) number of bins}
  \item{names}{variable names}   
}

\details{
  If \code{binw} is missing, the default binwidth is \eqn{\hat{b}_i = 2 \cdot 3^{1/(d+2)} \pi^{d/(2d+4)} S_i
  n^{-1/(d+2)}}{b_i = 2*3^(1/(d+2))*pi^(d/(2d+4))*S_i*n^(-1/(d+2))}, the
  normal scale selector.  

  If \code{xmin} is missing then it defaults to the data minimum.  If
  \code{xmax} is missing then it defaults to the data maximum.  
}

\seealso{\code{\link{plot.histde}}}

\examples{
## positive data example
set.seed(8192)
x <- 2^rnorm(100)
fhat <- histde(x=x)
plot(fhat, border=6)
points(c(0.5, 1), predict(fhat, x=c(0.5, 1)))

## large data example on a non-default grid
set.seed(8192)
x <- rmvnorm.mixt(10000, mus=c(0,0), Sigmas=invvech(c(1,0.8,1)))
fhat <- histde(x=x, xmin=c(-5,-5), xmax=c(5,5))
plot(fhat)

## See other examples in ? plot.histde
}

\keyword{smooth}
