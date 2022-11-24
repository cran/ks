\name{vector}
\alias{vec}
\alias{vech}
\alias{invvec}
\alias{invvech}
\title{Vector and vector half operators}

\description{
The vec (vector) operator takes a \eqn{d \times d}{d x d} matrix and stacks the
columns into a single vector of length \eqn{d^2}{d^2}. The vech (vector
half) operator
takes a symmetric \eqn{d \times d}{d x d} matrix and stacks the lower
triangular half into a single vector of length \eqn{d(d+1)/2}{d(d+1)/2}.
The functions invvec and invvech are the inverses of vec and
vech i.e. they form matrices from vectors. 
}

\usage{
vec(x, byrow=FALSE)
vech(x)
invvec(x, ncol, nrow, byrow=FALSE)
invvech(x)
}

\arguments{
  \item{x}{vector or matrix}
  \item{ncol,nrow}{number of columns and rows for inverse of vech}
  \item{byrow}{flag for stacking row-wise or column-wise. Default is FALSE.}  
  }
  
\references{ Magnus, J.R. & Neudecker H.M. (2007) \emph{Matrix
  Differential Calculus with Applications in Statistics and
  Econometrics (3rd edition)}, Wiley & Sons. Chichester.
}

\examples{
x <- matrix(1:9, nrow=3, ncol=3)
vec(x)
invvec(vec(x))
}

\keyword{algebra}
