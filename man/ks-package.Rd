\name{ks}
\alias{ks}
\docType{package}
\title{ks}
\description{
  Kernel smoothing for data from 1- to 6-dimensions.
}
\details{
  There are three main types of functions in this package:
  (a) computing kernel estimators, (b) computing bandwidth selectors and (c) displaying kernel estimators. The kernel used throughout is the normal (Gaussian) kernel.

  --For kernel density estimation, \code{\link{kde}} computes  
  \deqn{\hat{f}(\bold{x}; \bold{{\rm H}}) = n^{-1} \sum_{i=1}^n K_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i).}{hat(f)(x;H) = n^(-1) sum_i K_H (x - X_i).} 

  For display, its \code{plot} method calls \code{\link{plot.kde}}.

  --For kernel discriminant analysis, \code{\link{kda.kde}} computes density estimates for each the 
  groups in the training data, and the discriminant surface. 
  Its \code{plot} method is \code{\link{plot.kda.kde}}. The wrapper function \code{\link{Hkda}} computes 
  bandwidth matrices for each group in the training data, by calling the above selectors.
  
  --For kernel density derivative estimation, the main function is \code{\link{kdde}}
  \deqn{\widehat{{\sf D}^{\otimes r}f}(\bold{x}; \bold{{\rm H}}) = n^{-1} \sum_{i=1}^n {\sf D}^{\otimes r}K_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i).}{hat(f)^(r)(x;H) = n^(-1) sum_i D^r K_H (x - X_i).} 

  --There are several varieties of bandwidth matrix selectors 
  \itemize{
  \item plug-in \code{\link{hpi}} (1-d); 
  \code{\link{Hpi}}, \code{\link{Hpi.diag}} (2- to 6-d) \cr
  \item least squares (or unbiased) cross validation (LSCV or UCV) \code{\link{hlscv}} (1-d);
  \code{\link{Hlscv}}, \code{\link{Hlscv.diag}} (2- to 6-d) \cr
  \item biased cross validation (BCV) 
  \code{\link{Hbcv}}, \code{\link{Hbcv.diag}} (2- to 6-d) \cr
  \item smoothed cross validation (SCV) \code{\link{hscv}} (1-d);
  \code{\link{Hscv}}, \code{\link{Hscv.diag}} (2- to 6-d) \cr
  \item normal scale \code{\link{hmise.mixt}},
  \code{\link{hamise.mixt}} (1-d); and \code{\link{Hmise.mixt}},
  \code{\link{Hamise.mixt}}, \code{\link{Hmise.mixt.diag}}, \code{\link{Hamise.mixt.diag}} (2- to 6-d).}
  For 1-d data, the bandwidth \eqn{h}{h} is the standard deviation of the normal kernel, whereas for multivariate data, the bandwidth matrix 
  \eqn{\bold{{\rm H}}}{H} is the variance matrix.

  --For kernel functional estimation, \code{kfe} computes the \eqn{r}{r}-th order integrated density functional  
  \deqn{\hat{{\bold \psi}}_r (\bold{{\rm H}}) = n^{-2} \sum_{i=1}^n \sum_{j=1}^n {\sf D}^{\otimes r}K_{\bold{{\rm H}}}(\bold{X}_i-\bold{X}_j).}{hat(psi)_r = n^(-2) sum_i sum_j D^r K_H (X_i - X_j).} The plug-in selector is \code{\link{Hpi.kfe}}.

  --Binned kernel estimation is available for d = 1, 2, 3, 4. This makes kernel estimators 
  feasible for large samples.
  
  --For an overview of this package with 2-d density estimation, see 
  \code{vignette("kde")}. 
}
\author{
  Tarn Duong for most of the package. 
  M.P. Wand for the binned estimation, univariate plug-in selector
  and univariate density derivative estimator code.
  Jose E. Chacon for the unconstrained pilot
  functional estimation and fast implementation of density derivative estimation code.   
}

\references{
  Bowman, A. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Oxford University Press, Oxford.
  
  Duong, T. (2004) \emph{Bandwidth Matrices for Multivariate Kernel Density 
    Estimation.} Ph.D. Thesis, University of Western Australia. 

  Scott, D.W. (1992) \emph{Multivariate Density Estimation: Theory,
    Practice, and Visualization}. John Wiley & Sons, New York.

  Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman & Hall/CRC, London.

  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag. New York.

  Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman &
  Hall/CRC, London.
}

\keyword{package}

\seealso{\code{sm}, \code{KernSmooth}
}

