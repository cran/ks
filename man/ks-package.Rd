\name{ks}
\alias{ks}
\docType{package}
\title{
ks
}
\description{
  Kernel density estimation and kernel discriminant analysis for
  data from 1- to 6-dimensions, with display functions.
  
}
\details{
  There are three main types of functions in this package:
  (a) computing bandwidth
  selectors, (b) computing kernel estimators and (c) displaying kernel estimators.

  \enumerate{
	\item[(a)] For the bandwidth matrix selectors, there are several varieties:
  \enumerate{
  \item[(i)] plug-in \code{\link{hpi}} (1-d); 
  \code{\link{Hpi}}, \code{\link{Hpi.diag}} (2- to 6-d)
  \item[(ii)] least squares (or unbiased) cross validation (LSCV or UCV)
  \code{\link{Hlscv}}, \code{\link{Hlscv.diag}} (2- to 6-d)
  \item[(iii)] biased cross validation (BCV) 
  \code{\link{Hbcv}}, \code{\link{Hbcv.diag}} (2- to 6-d)
  \item[(iv)] smoothed cross validation (SCV)
  \code{\link{Hscv}}, \code{\link{Hscv.diag}} (2- to 6-d)
  \item[(v)] normal scale selectors \code{\link{hmise.mixt}},
  \code{\link{hamise.mixt}} (1-d); and \code{\link{Hmise.mixt}},
  \code{\link{Hamise.mixt}} (2- to 6-d).
  }
  For multivariate selectors, the first selector in each pair is the
  unconstrained (or full) 
  selector, and the second one is the diagonal selectors.
  
  \item[(b)] For kernel density estimation, the main function is
  \code{\link{kde}}. For kernel discriminant analysis, it's \code{\link{kda.kde}}.

  \item[(c)] For display, \code{plot} via (\code{\link{plot.kde}} and
  \code{\link{plot.kda.kde}}) sends to a graphics window 
  the results of density estimation or discriminant analysis.
  }

  For d = 1, 2, 3, 4, binned kernel estimation is available. 
  
  For an overview of this package with 2-d density estimation, see 
  \code{vignette("kde")}. 
}
\author{
  Tarn Duong <tarn.duong@gmail.com> for most of the package. 
  Matt Wand for the binned estimation and univariate plug-in selector
  code. Jose E. Chac\'on for the unconstrained pilot
  functional estimation and (A)MISE-optimal selectors for normal mixture
  densities code.   
}
\references{
  Bowman, A. \& Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Oxford University Press. Oxford.

  Chac\'on, J.E. \& Duong, T. (2008) Multivariate plug-in bandwidth
  selection with unconstrained pilot matrices. \emph{Submitted.}
  
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

  Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman \& Hall/CRC. London.

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

