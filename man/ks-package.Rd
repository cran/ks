\name{ks-package}
\alias{ks}
\docType{package}
\title{ks}
\description{
  Kernel smoothing for data from 1- to 6-dimensions.
}

\details{
  There are three main types of functions in this package:
  \itemize{
    \item computing kernel estimators - these function names begin with `k'
    \item computing bandwidth selectors - these begin with `h' (1-d) or
    `H' (>1-d)
    \item displaying kernel estimators - these begin with `plot'.
  }
  The kernel used throughout is the normal (Gaussian) kernel \eqn{K}{K}.
  For 1-d data, the bandwidth \eqn{h}{h} is the standard deviation of
  the normal kernel, whereas for multivariate data, the bandwidth matrix 
  \eqn{\bold{{\rm H}}}{H} is the variance matrix.
 
  --For kernel density estimation, \code{\link{kde}} computes  
  \deqn{\hat{f}(\bold{x}) = n^{-1} \sum_{i=1}^n K_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i).}{hat(f)(x) = n^(-1) sum_i K_H (x - X_i).} 

  The bandwidth matrix \eqn{\bold{{\rm H}}}{H} is a matrix of smoothing
  parameters and its choice is crucial for the performance of kernel
  estimators. For display, its \code{plot} method calls \code{\link{plot.kde}}.

  --For kernel density estimators, there are several varieties of bandwidth selectors 
  \itemize{
  \item plug-in \code{\link{hpi}} (1-d); 
  \code{\link{Hpi}}, \code{\link{Hpi.diag}} (2- to 6-d) 
  \item least squares (or unbiased) cross validation (LSCV or UCV) \code{\link{hlscv}} (1-d);
  \code{\link{Hlscv}}, \code{\link{Hlscv.diag}} (2- to 6-d) 
  \item biased cross validation (BCV) 
  \code{\link{Hbcv}}, \code{\link{Hbcv.diag}} (2- to 6-d) 
  \item smoothed cross validation (SCV) \code{\link{hscv}} (1-d);
  \code{\link{Hscv}}, \code{\link{Hscv.diag}} (2- to 6-d) 
  \item normal scale \code{\link{hns}} (1-d); \code{\link{Hns}} (2- to 6-d).
  }
  
  --For kernel density derivative estimation, the main function is \code{\link{kdde}}
  \deqn{{\sf D}^{\otimes r}\hat{f}(\bold{x}) = n^{-1} \sum_{i=1}^n {\sf
  D}^{\otimes r}K_{\bold{{\rm H}}} (\bold{x} -
\bold{X}_i).}{hat(f)^(r)(x) = n^(-1) sum_i D^r K_H (x - X_i).} 
  The bandwidth selectors are a modified subset of those for
  \code{\link{kde}}, i.e. \code{\link{Hlscv}}, \code{\link{Hns}}, \code{\link{Hpi}}, \code{\link{Hscv}} with \code{deriv.order>0}. 
   Its \code{plot} method is \code{\link{plot.kdde}} for plotting each
   partial derivative singly.

  --For kernel discriminant analysis, the main function is
  \code{\link{kda}} which computes density estimates for each the 
  groups in the training data, and the discriminant surface. 
  Its \code{plot} method is \code{\link{plot.kda}}. The wrapper function
  \code{\link{hkda}}, \code{\link{Hkda}} computes 
  bandwidths for each group in the training data for \code{kde},
  e.g. \code{hpi}, \code{Hpi}.
  
  --For kernel functional estimation, the main function is
  \code{kfe} which computes the \eqn{r}{r}-th order integrated density functional  
  \deqn{\hat{{\bold \psi}}_r  = n^{-2} \sum_{i=1}^n \sum_{j=1}^n {\sf D}^{\otimes r}K_{\bold{{\rm H}}}(\bold{X}_i-\bold{X}_j).}{hat(psi)_r = n^(-2) sum_i sum_j D^r K_H (X_i - X_j).} The plug-in selectors are \code{\link{hpi.kfe}} (1-d), \code{\link{Hpi.kfe}} (2- to 6-d).
  Kernel functional estimates are usually not required to computed
  directly by the user, but only within other functions in the package.
    
  --For kernel-based 2-sample testing, the main function is
  \code{\link{kde.test}} which computes the integrated 
  \eqn{L_2}{L2} distance between the two density estimates as the test
  statistic, comprising a linear combination of 0-th order kernel
  functional estimates:
  \deqn{\hat{T} = \hat{\psi}_{0,1} + \hat{\psi}_{0,2} - (\hat{\psi}_{0,12} +
  \hat{\psi}_{0,21}),}{hat(T) = hat(psi)_0,1 + hat(psi)_0,2 - (hat(psi)_0,12 +
  hat(psi)_0,21),} and the corresponding p-value. The \eqn{\psi}{psi} are
  zero order kernel functional estimates with the subscripts indicating
  that 1 = sample 1 only, 2 = sample 2 only, and 12, 21 =
  samples 1 and 2.  The bandwidth selectors are \code{\link{hpi.kfe}},
  \code{\link{Hpi.kfe}} with \code{deriv.order=0}. 

  --For kernel-based local 2-sample testing, the main function is
   \code{\link{kde.local.test}} which computes the squared distance
   between the two density estimates as the test  
   statistic \deqn{\hat{U}(\bold{x}) = [\hat{f}_1(\bold{x}) - \hat{f}_2(\bold{x})]^2}{hat(U)(x) =
   [hat(f)_1(x) - hat(f)_2(x)]^2} and the corresponding local
   p-values.  The bandwidth selectors are those used with \code{\link{kde}},
   e.g. \code{\link{hpi}, \link{Hpi}}.  

  --For kernel cumulative distribution function estimation, the main
  function is \code{\link{kcde}}
  \deqn{\hat{F}(\bold{x}) = n^{-1} \sum_{i=1}^n
  \mathcal{K}_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i)}{hat(F)(x) =
  n^(-1) sum_i intK_H (x - X_i)}
  where \eqn{\mathcal{K}}{intK} is the integrated kernel.  
  The bandwidth selectors are \code{\link{hpi.kcde}},
  \code{\link{Hpi.kcde}}.  Its \code{plot} method is
  \code{\link{plot.kcde}}.
  There exist analogous functions for the survival function \eqn{\hat{\bar{F}}}{hat(bar(F))}.

  --For kernel estimation of a ROC (receiver operating characteristic)
  curve to compare two samples from \eqn{\hat{F}_1,
  \hat{F}_2}{hat(F)_1, hat(F)_2}, the main function is \code{\link{kroc}}
  \deqn{(\hat{F}_{\hat{Y}_1}(z),
  \hat{F}_{\hat{Y}_2}(z))}{(hat(F)_hat(Y1))(z), hat(F_hat(Y2))(z))} based on the cumulative distribution functions of
  \eqn{\hat{Y}_j = \hat{\bar{F}}_1(\bold{X}_j), j=1,2}{hat(Yj)=hat(bar(F))_1(X_j), j=1,2}.
    
  The bandwidth selectors are those used with \code{\link{kcde}},
  e.g. \code{\link{hpi.kcde}, \link{Hpi.kcde}}  for
  \eqn{\hat{F}_{\hat{Y}_j}, \hat{\bar{F}}_1}{hat(F)_hat(Yj), hat(bar(F))_1}. Its \code{plot} method
  is \code{\link{plot.kroc}}.   

  --For kernel estimation of a copula, the
  main function is \code{\link{kcopula}}
  \deqn{\hat{C}(\bold{z}) = \hat{F}(\hat{F}_1^{-1}(z_1), \dots,
  \hat{F}_d^{-1}(z_d))}{hat(C)(z) = hat(F)(hat(F)_1^(-1)(z_1),..., hat(F)_d^(-1)(z_d))}
  where \eqn{\hat{F}_j^{-1}(z_j)}{hat(F)_j^(-1)(z_j)} is
  the \eqn{z_j}{z_j}-th quantile of of the \eqn{j}{j}-th marginal
  distribution \eqn{\hat{F}_j}{hat(F_j)}.
  The bandwidth selectors are those used with \code{\link{kcde}} for
  \eqn{\hat{F}, \hat{F}_j}{hat(F), hat(F)_j}.
  Its \code{plot} method is \code{\link{plot.kcde}}.

  --For kernel estimation of a copula density, the
  main function is \code{\link{kcopula.de}}
  \deqn{\hat{c}(\bold{z}) = n^{-1} \sum_{i=1}^n
  K_{\bold{{\rm H}}} (\bold{z} - \hat{\bold{Z}}_i)}{hat(c)(z) =
  hat(f)(z) = n^(-1) sum_i K_H (z - hat(Z)_i)}
  where \eqn{\hat{\bold{Z}}_i = (\hat{F}_1(X_{i1}), \dots,
    \hat{F}_d(X_{id}))}{hat(Z)_i = (hat(F)_1(X_i1), \dots, hat(F)_d(X_id))}.
  The bandwidth selectors are those used with \code{\link{kde}} for
  \eqn{\hat{c}}{hat(c)} and \code{\link{kcde}} for \eqn{\hat{F}_j}{hat(F)_j}.
  Its \code{plot} method is \code{\link{plot.kde}}.
  
  --Binned kernel estimation is available for d = 1, 2, 3, 4. This makes kernel estimators feasible for large samples.
  
  --For an overview of this package with 2-d density estimation, see 
  \code{vignette("kde")}. 
}
\author{
  Tarn Duong for most of the package. 
  M.P. Wand for the binned estimation, univariate plug-in selector
  and univariate density derivative estimator code.
  Jose E. Chacon for the unconstrained pilot
  functional estimation and fast implementation of derivative-based
  estimation code.
  Artur and Jaroslaw Gramacki for the binned estimation for
  unconstrained bandwidth matrices. 
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

