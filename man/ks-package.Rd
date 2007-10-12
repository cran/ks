\name{ks}
\alias{ks}
\docType{package}
\title{
ks
}
\description{
  Kernel density estimation and kernel discriminant analysis for
  multivariate data (1- to 6-dimensions) with display functions.
  
}
\details{
  There are three main types of functions in this package: (a) computing bandwidth
  selectors, (b) computing kernel estimators and (c) displaying kernel estimators.

  (a) For the bandwidth matrix selectors, there are several varieties:
  (i) plug-in
  \code{\link{Hpi}}, \code{\link{Hpi.diag}}
  (ii) least squares (or unbiased) cross validation (LSCV or UCV)
  \code{\link{Hlscv}}, \code{\link{Hlscv.diag}},
  (iii) biased cross validation (BCV)
  \code{\link{Hbcv}}, \code{\link{Hbcv.diag}} and
  (iv) smoothed cross validation (SCV)
  \code{\link{Hscv}}, \code{\link{Hscv.diag}}.
	The first selector in each pair is the unconstrained (or full)
	selector,
	and the sceond one is the diagonal selectors.

	Scalar bandwidth selectors are not provided - see
  \code{sm} or \code{KernSmooth} packages. 

  (b) For kernel density estimation, the main function is
  \code{\link{kde}}. For kernel discriminant analysis, it's \code{\link{kda.kde}}.

  (c) For display, versions of \code{plot}, \code{\link{plot.kde}} and
  \code{\link{plot.kda.kde}}, send to a graphics window 
  the results of density estimation or discriminant analysis.

  For d = 1, 2, 3, 4, binned kernel estimation is available. 
  
  For an overview of this package with 2-dim density estimation, see 
  \code{vignette("kde")}. 
}
\author{
Tarn Duong <tduong@pasteur.fr> for most of the package. 
Matt Wand for the binned estimation code. Jose E. Chacon for the pilot
functional estimation code.  
}
\references{
  Bowman, A. \& Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Oxford University Press. Oxford.
  
  Duong, T. (2004) \emph{Bandwidth Matrices for Multivariate Kernel Density 
    Estimation.} Ph.D. Thesis. University of Western Australia. 
 
  Duong, T. \& Hazelton, M.L. (2003) Plug-in bandwidth matrices for
    bivariate kernel density estimation. \emph{Journal of Nonparametric
  Statistics}, \bold{15}, 17-30.
  
  Duong, T. \& Hazelton, M.L. (2005) Cross-validation bandwidth
    matrices for multivariate kernel density estimation. \emph{Scandinavian Journal 
  of Statistics}, \bold{32}, 485-506.

  Sain, S.R., Baggerly, K.A. \& Scott, D.W. (1994)
  Cross-validation of multivariate densities. \emph{Journal of the
  American Statistical Association}. \bold{82}, 1131-1146.

  Scott, D.W. (1992) \emph{Multivariate Density Estimation: Theory,
    Practice, and Visualization}. John Wiley \& Sons. New York.
  
  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.

  Wand, M.P. \& Jones, M.C. (1994) Multivariate plugin bandwidth
    selection. \emph{Computational Statistics}, \bold{9}, 97-116.
  
  Wand, M.P. \& Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman \&
  Hall/CRC. London.
}

\keyword{ package }

\seealso{
\code{\link[sm:sm-package]{sm}}, \code{\link[KernSmooth]{KernSmooth}}
}

