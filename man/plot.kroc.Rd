\name{plot.kroc}
\alias{plot.kroc}
\title{Plot for kernel receiver operating characteristic curve (ROC) estimate}

\description{
  Plot for kernel receiver operating characteristic curve (ROC) estimate 1- to 3-dimensional data.
}

\usage{
\method{plot}{kroc}(x, add=FALSE, add.roc.ref=FALSE, xlab, ylab, 
   alpha=1, col=1, ...)
}

\arguments{
  \item{x}{object of class \code{kroc} (output from \code{\link{kroc}})}
  \item{add}{flag to add to current plot. Default is FALSE.}
  \item{add.roc.ref}{flag to add reference ROC curve. Default is FALSE.} 
  \item{xlab}{x-axis label. Default is "False positive rate
  (bar(specificity))".}
  \item{ylab}{y-axis label. Default is "True positive rate (sensitivity)".} 
  \item{alpha,col}{transparency value and colour of line}
  \item{...}{other graphics parameters used in \code{\link{plot.kde}}.} 
}
  
\value{
  Plots for 1-d and 2-d are sent to graphics window. Plot for 3-d is sent to
  graphics/RGL window. 
}

\seealso{\code{\link{plot.kde}}}

\examples{
library(MASS)
data(fgl)
x1 <- fgl[fgl[,"type"]=="WinF",c("RI", "Na")]
x2 <- fgl[fgl[,"type"]=="Head",c("RI", "Na")]
Rhat <- kroc(x1=x1, x2=x2) 
plot(Rhat, add.roc.ref=TRUE)
}

\keyword{hplot}
