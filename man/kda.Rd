\name{kda, pda, compare}
\alias{kda}
\alias{pda}
\alias{compare}

\title{Kernel and parametric discriminant analysis}
\description{
  Kernel and parametric discriminant analysis.
}
\usage{
kda(x, x.group, Hs, y, prior.prob=NULL)
pda(x, x.group, y, prior.prob=NULL, type="quad")
compare(x.group, est.group)
}

\arguments{
  \item{x}{matrix of training data values}
  \item{x.group}{vector of group labels for training data}
  \item{est.group}{vector of estimated group labels}
  \item{y}{matrix of test data} 
  \item{Hs}{(stacked) matrix of bandwidth matrices}
  \item{prior.prob}{vector of prior probabilities}
  \item{type}{\code{"line"} = linear discriminant, \code{"quad"} =
    quadratic discriminant}
}

\value{
  The discriminant analysers are \code{kda} and \code{pda} and these
  return a vector of group labels assigned via discriminant
  analysis.  If the test data \code{y} are given then these are
  classified. Otherwise the training data \code{x} are classified.
  
  The function \code{compare} creates a comparison between the true
  group labels \code{x.group} and the estimated ones \code{est.group}.
  It returns a list with fields
  \item{cross}{cross-classification table with the rows
    indicating the true group and the columns the estimated group}
  \item{error}{misclassification rate (MR) where 
    \deqn{\textrm{MR} = \frac{\textrm{number of points wrongly classified}}{\textrm{total number of
	points}}}{MR = (number of points wrongly classified) / (total number of
    points)}}
Note that this MR is only suitable when we have test data. If  we
don't have test data, then the cross validated estimate is more
appropriate.  See Silverman (1986).
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
  \code{MASS} library namely \code{lda} and
  \code{qda}
  for linear and quadratic discriminants.
}

\seealso{  \code{\link{kda.kde}}, \code{\link{pda.pde}}}
}

\examples{

### bivariate example - restricted iris dataset  
library(MASS)
data(iris)
iris.mat <- rbind(iris[,,1], iris[,,2], iris[,,3])
ir <- iris.mat[,c(1,2)]
ir.gr <- iris.mat[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, H, ir)
lda.gr <- pda(ir, ir.gr, ir, type="line")
qda.gr <- pda(ir, ir.gr, ir, type="quad")

compare(kda.gr, ir.gr)
compare(qda.gr, ir.gr)
compare(lda.gr, ir.gr)

### multivariate example - full iris dataset
ir <- iris[,1:4]
ir.gr <- iris[,5]

H <- Hkda(ir, ir.gr, bw="plugin", pre="scale")
kda.gr <- kda(ir, ir.gr, H, ir)
lda.gr <- pda(ir, ir.gr, ir, type="line")
qda.gr <- pda(ir, ir.gr, ir, type="quad")

compare(kda.gr, ir.gr)
compare(qda.gr, ir.gr)
compare(lda.gr, ir.gr)
}

\keyword{ smooth }
