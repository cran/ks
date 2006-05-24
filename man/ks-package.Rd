\name{ks}
\alias{ks}
\docType{package}
\title{
ks
}
\description{
  Kernel density estimators and kernel discriminant analysis for
  multivariate data (2- to 6-dimensions) with display functions.
  
}
\details{
  There are three main types of functions in this package: (a) bandwidth
  selectors, (b) kernel density estimators and (c) display.

  For the bandwidth matrix selectors, there are several varieties: (a) plug-in
  \code{\link{Hpi}}, (b) least squares (or unbiased) cross validation
  (LSCV or UCV)
  \code{\link{Hlscv}}, (c) biased cross validation (BCV)
  \code{\link{Hbcv}} and (d) smoothed cross validation (SCV)
  \code{\link{Hscv}}.

  For kernel density estimation, the main function is
  \code{\link{kde}}. For kernel discriminant analysis, it's \code{\link{kda}}.

  For display, versions of \code{plot} send to a graphics window (the
  usual one for 2-dimensions or an RGL window for 3-dimensions)
  the results of density estimation or discriminant analysis.

  For 1- to 4-dimensional data, binned versions of the kernel density
  estimators are also available, adapted from the
  \code{\link[KernSmooth]{KernSmooth}} library. These binned estimators are 
  useful for large sample sizes. 
  
}
\author{
Tarn Duong <tduong@maths.unsw.edu.au>
}
\references{
  Bowman, A. \& Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Oxford University Press. Oxford.
  
  Duong, T. (2004) \emph{Bandwidth Matrices for Multivariate Kernel Density 
    Estimation.} Ph.D. Thesis. University of Western Australia. 
 
  Duong, T. \& Hazelton, M.L. (2003) \emph{Plug-in bandwidth matrices for
    bivariate kernel density estimation}. Journal of Nonparametric
  Statistics \bold{15}, 17-30.
  
  Duong, T. \& Hazelton, M.L. (2005) \emph{Cross-validation bandwidth
    matrices for multivariate kernel density estimation}. Scandinavian Journal 
  of Statistics. \bold{32}, 485-506.

  Sain, S.R., Baggerly, K.A. \& Scott, D.W. (1994)
  \emph{Cross-validation of multivariate densities}. Journal of the
  American Statistical Association. \bold{82}, 1131-1146.

  Scott, D.W. (1992) \emph{Multivariate Density Estimation: Theory,
    Practice, and Visualization}. John Wiley \& Sons. Mew York.
  
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.

  Wand, M.P. \& Jones, M.C. (1994) \emph{Multivariate plugin bandwidth
    selection}.  Computational Statistics \bold{9}, 97-116.
  
  Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman \&
  Hall/CRC. London.
}

\keyword{ package }

\seealso{
\code{\link[sm:sm-package]{sm}}, \code{\link[KernSmooth]{KernSmooth}}
}

