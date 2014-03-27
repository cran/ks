\name{plot.kroc}
\alias{plot.kroc}

\title{Plot for kernel receiver operating
  characteristic curve (ROC) estimate}
\description{
  Plot for kernel receiver operating
  characteristic curve (ROC) estimate 1- to 3-dimensional data.
}

\usage{
\method{plot}{kroc}(x, add=FALSE, add.roc.ref=FALSE, ylab, xlab, ...)
}

\arguments{
  \item{x}{an object of class \code{kroc} (output from \code{\link{kroc}})}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{add.roc.ref}{flag to add reference ROC curve. Default is FALSE.} 
  \item{xlab}{x-axis label. Default is "False positive rate
  (bar(specificity))".}
  \item{ylab}{y-axis label. Default is "True positive rate (sensitivity)".} 
  \item{...}{other graphics parameters used in \code{\link{plot.kde}}.} 
}
  

\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  RGL window. 
}


\seealso{\code{\link{plot.kde}}}

\examples{
library(MASS)
data(iris)
x1 <- iris[iris[,5]=="versicolor",1:2]
x2 <- iris[iris[,5]=="virginica",1:2]
Rhat <- kroc(x1=x1, x2=x2, nref=1e3) 
plot(Rhat)
abline(0,1, lty=2, col="grey")
}


\keyword{hplot}
