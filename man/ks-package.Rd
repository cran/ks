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
    
  --For kernel density estimation, there are several varieties of bandwidth selectors 
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

  --For kernel density support estimation, the main function is
  \code{\link{ksupp}} which is (the convex hull of) 
  \deqn{\{\bold{x}: \hat{f}(\bold{x}) >
    \tau\}}{\{x: hat(f) > tau\}}
  for a suitable level \eqn{\tau}{tau}. This is closely related to the \eqn{\tau}{tau}-level set of
  \eqn{\hat{f}}{hat(f)}. 
  
  --For truncated kernel density estimation, the main function is
  \code{\link{kde.truncate}}
  \deqn{\hat{f} (\bold{x}) \bold{1}\{\bold{x} \in \Omega\} /
  \int_{\Omega}\hat{f} (\bold{x}) \, d\bold{x}}{hat(f)(x) 1\{x in Omega\}/int hat(f) 1\{x in Omega\}}
  for a bounded data support \eqn{\Omega}{Omega}. The standard density
  estimate \eqn{\hat{f}}{hat(f)} is truncated and rescaled to give
  unit integral over \eqn{\Omega}{Omega}. Its \code{plot} method calls \code{\link{plot.kde}}.

  --For boundary kernel density estimation where the kernel function is
  modified explicitly in the boundary region, the main function is
  \code{\link{kde.boundary}}
  \deqn{ n^{-1} \sum_{i=1}^n K^*_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i)}{hat(f)(x) = n^(-1) sum_i K*_H (x - X_i)} 
  for a boundary kernel \eqn{K^*}{K*}. Its \code{plot} method calls \code{\link{plot.kde}}.

  --For variable kernel density estimation where the bandwidth is not a
  constant matrix, the main functions are \code{\link{kde.balloon}}
  \deqn{\hat{f}_{\rm ball}(\bold{x}) = n^{-1} \sum_{i=1}^n K_{\bold{{\rm H}}(\bold{x})} (\bold{x} - \bold{X}_i)}{hat(f)_ball(x) = n^(-1) sum_i K_H(x) (x - X_i)}
  and
  \code{\link{kde.sp}} 
  \deqn{\hat{f}_{\rm SP}(\bold{x}) = n^{-1} \sum_{i=1}^n K_{\bold{{\rm H}}(\bold{X}_i)} (\bold{x} - \bold{X}_i).}{hat(f)_SP(x) = n^(-1) sum_i K_H(X_i) (x - X_i).}
  For the balloon estimation \eqn{\hat{f}_{\rm ball}}{hat(f)_ball} the
  bandwidth varies with the estimation point \eqn{\bold{x}}{x}, whereas
  for the sample point estimation  \eqn{\hat{f}_{\rm SP}}{hat(f)_SP}
  the bandwidth varies with the data point
  \eqn{\bold{X}_i, i=1,\dots,n}{X_i, i=1, ..., n}. 
  Their \code{plot} methods call \code{\link{plot.kde}}. The bandwidth
  selectors for \code{kde.balloon} are based on the normal scale bandwidth
  \code{Hns(,deriv.order=2)} via the MSE minimal formula, and for
  \code{kde.SP} on \code{Hns(,deriv.order=4)} via the Abramson formula.
  
  --For kernel density derivative estimation, the main function is \code{\link{kdde}}
  \deqn{{\sf D}^{\otimes r}\hat{f}(\bold{x}) = n^{-1} \sum_{i=1}^n {\sf
  D}^{\otimes r}K_{\bold{{\rm H}}} (\bold{x} -
\bold{X}_i).}{hat(f)^(r)(x) = n^(-1) sum_i D^r K_H (x - X_i).} 
  The bandwidth selectors are a modified subset of those for
  \code{\link{kde}}, i.e. \code{\link{Hlscv}}, \code{\link{Hns}}, \code{\link{Hpi}}, \code{\link{Hscv}} with \code{deriv.order>0}. 
   Its \code{plot} method is \code{\link{plot.kdde}} for plotting each
   partial derivative singly.

  --For kernel summary curvature estimation, the main function is
  \code{\link{kcurv}}
  \deqn{\hat{s}(\bold{x})= - \bold{1}\{{\sf D}^2 \hat{f}(\bold{x}) <
    0\} \mathrm{abs}(|{\sf D}^2 \hat{f}(\bold{x})|)}{hat(s)(x) = -1\{D^2 hat(f)(x) <0)*abs(det(D^2 hat(f)(x)))\}} where \eqn{{\sf D}^2
    \hat{f}(\bold{x})}{D^2 hat(f)(x)} is the kernel Hessian matrix estimate.
  It has the same structure as a kernel density estimate so its \code{plot}
  method calls \code{\link{plot.kde}}. 
   
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
   statistic \deqn{\hat{U}(\bold{x}) = [\hat{f}_1(\bold{x}) -
     \hat{f}_2(\bold{x})]^2}{hat(U)(x) = [hat(f)_1(x) - hat(f)_2(x)]^2} and the corresponding local
   p-values.  The bandwidth selectors are those used with \code{\link{kde}},
   e.g. \code{\link{hpi}, \link{Hpi}}.  

  --For kernel cumulative distribution function estimation, the main
  function is \code{\link{kcde}}
  \deqn{\hat{F}(\bold{x}) = n^{-1} \sum_{i=1}^n
    \mathcal{K}_{\bold{{\rm H}}} (\bold{x} - \bold{X}_i)}{hat(F)(x) = n^(-1) sum_i intK_H (x - X_i)}
  where \eqn{\mathcal{K}}{intK} is the integrated kernel.  
  The bandwidth selectors are \code{\link{hpi.kcde}},
  \code{\link{Hpi.kcde}}.  Its \code{plot} method is
  \code{\link{plot.kcde}}.
  There exist analogous functions for the survival function \eqn{\hat{\bar{F}}}{hat(bar(F))}.

  --For kernel estimation of a ROC (receiver operating characteristic)
  curve to compare two samples from \eqn{\hat{F}_1,
  \hat{F}_2}{hat(F)_1, hat(F)_2}, the main function is \code{\link{kroc}}
  \deqn{\{\hat{F}_{\hat{Y}_1}(z),
  \hat{F}_{\hat{Y}_2}(z)\}}{\{hat(F)_hat(Y1))(z), hat(F_hat(Y2))(z)\}} based on the cumulative distribution functions of
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

  %--For kernel estimation of a copula density, the
  %main function is \code{\link{kcopula.de}}
  %\deqn{\hat{c}(\bold{z}) = n^{-1} \sum_{i=1}^n
  %K_{\bold{{\rm H}}} (\bold{z} - \hat{\bold{Z}}_i)}{hat(c)(z) =
  %hat(f)(z) = n^(-1) sum_i K_H (z - hat(Z)_i)}
  %where \eqn{\hat{\bold{Z}}_i = (\hat{F}_1(X_{i1}), \dots,
  %  \hat{F}_d(X_{id}))}{hat(Z)_i = (hat(F)_1(X_i1), \dots, hat(F)_d(X_id))}.
  %The bandwidth selectors are those used with \code{\link{kde}} for
  %\eqn{\hat{c}}{hat(c)} and \code{\link{kcde}} for \eqn{\hat{F}_j}{hat(F)_j}.
  %Its \code{plot} method is \code{\link{plot.kde}}.

  --For kernel mean shift clustering, the main function is
  \code{\link{kms}}. The mean shift recurrence relation of the candidate
  point \eqn{{\bold x}}{x}   
  \deqn{{\bold x}_{j+1} = {\bold x}_j + \bold{{\rm H}} {\sf D} \hat{f}({\bold
      x}_j)/\hat{f}({\bold x}_j),}{x_j+1 = x_j + H D hat(f)(x_j)/hat(f)(x_j),}
  where \eqn{j\geq 0}{j>=0} and \eqn{{\bold x}_0 = {\bold x}}{x_0 = x},
  is iterated until \eqn{{\bold x}}{x} converges to its
  local mode in the density estimate \eqn{\hat{f}}{hat(f)} by following
  the density gradient ascent paths. This mode determines the cluster
  label for \eqn{\bold{x}}{x}.  The bandwidth selectors are those used with
  \code{\link{kdde}(,deriv.order=1)}.  

  --For kernel density ridge estimation, the main function is
    \code{\link{kdr}}. The kernel density ridge recurrence relation of
    the candidate point \eqn{{\bold x}}{x} 
    \deqn{{\bold x}_{j+1} = {\bold x}_j + \bold{{\rm U}}_{(d-1)}({\bold
    x}_j)\bold{{\rm U}}_{(d-1)}({\bold x}_j)^T \bold{{\rm H}} {\sf D}
    \hat{f}({\bold x}_j)/\hat{f}({\bold x}_j),}{x_j+1 = x_j +
    U_(d-1)(x_j) U_(d-1)(x_j)^T H D hat(f)(x_j)/hat(f)(x_j),}
    where \eqn{j\geq 0}{j>=0}, \eqn{{\bold x}_0 = {\bold x}}{x_0 =
    x} and \eqn{\bold{{\rm U}}_{(d-1)}}{U_(d-1)} is the 1-dimensional projected
    density gradient, 
    is iterated until \eqn{{\bold x}}{x} converges to the ridge in the
    density estimate. The bandwidth selectors are those used with 
    \code{\link{kdde}(,deriv.order=2)}. 
  
  -- For kernel feature significance, the main function
     \code{\link{kfs}}.  The hypothesis test at a point \eqn{\bold{x}}{x} is
  \eqn{H_0(\bold{x}): \mathsf{H} f(\bold{x}) < 0}{H0(x): H f(x) < 0},
     i.e. the density Hessian matrix \eqn{\mathsf{H} f(\bold{x})}{H f(x)} is negative definite.
  The test statistic is  
  \deqn{W(\bold{x}) = \Vert 
    \mathbf{S}(\bold{x})^{-1/2} \mathrm{vech} \ \mathsf{H} \hat{f} (\bold{x})\Vert ^2}{%
    W(x) = ||S(x)^(-1/2) vech H hat{f}(x)||^2}
  where \eqn{{\sf H}\hat{f}}{H hat{f}} is the
  Hessian estimate, vech is the vector-half operator, and
  \eqn{\mathbf{S}}{S} is an estimate of the null variance.
  \eqn{W(\bold{x})}{W(x)} is
  approximately \eqn{\chi^2}{chi-squared} distributed with
     \eqn{d(d+1)/2}{d(d+1)/2} degrees of freedom.
  If \eqn{H_0(\bold{x})}{H0(x)} is rejected, then \eqn{\bold{x}}{x}
  belongs to a significant modal region.
  The bandwidth selectors are those used with
  \code{\link{kdde}(,deriv.order=2)}. Its \code{plot} method is
  \code{\link{plot.kfs}}.
     
  --For deconvolution density estimation, the main function is
  \code{\link{kdcde}}. A weighted kernel density
  estimation with the contaminated data \eqn{{\bold W}_1, \dots, {\bold
  W}_n}{W_1, ..., W_n}, 
  \deqn{\hat{f}_w({\bold x}) = n^{-1} \sum_{i=1}^n
    \alpha_i K_{\bold{{\rm H}}}({\bold x} - {\bold W}_i),}{hat(f)(x) = n^(-1) sum_i alpha_i K_H (x - W_i),}
  is utilised, where the weights \eqn{\alpha_1, \dots,
  \alpha_n}{alpha_1, ..., alpha_n} are chosen via a
  quadratic optimisation involving the error variance and the
  regularisation parameter.  The bandwidth selectors are those used with
  \code{\link{kde}}.
     
  --Binned kernel estimation is an approximation to the exact kernel
  estimation and is available for d=1, 2, 3, 4. This makes
  kernel estimators feasible for large samples. 
  
  --For an overview of this package with 2-d density estimation, see 
  \code{vignette("kde")}.

  --For \pkg{ks} \eqn{\geq}{>=} 1.11.1, the \pkg{misc3d} and
  \pkg{rgl} (3-d plot), \pkg{OceanView} (quiver plot),
  \pkg{oz} (Australian map) packages have been moved from
  Depends to Suggests. This was done to allow \pkg{ks} to be installed
  on systems where these latter graphical-based packages can't be
  installed. Furthermore, since the future of OpenGL in R is not certain, 
  \pkg{plot3D} becomes the default for 3D plotting for \pkg{ks} \eqn{\geq}{>=} 1.12.0. 
  RGL plots are still supported though these may be deprecated in the 
  future.  
}

\author{
  Tarn Duong for most of the package. 
  M. P. Wand for the binned estimation, univariate plug-in selector
  and univariate density derivative estimator code.
  J. E. Chacon for the unconstrained pilot
  functional estimation and fast implementation of derivative-based
  estimation code.
  A. and J. Gramacki for the binned estimation for
  unconstrained bandwidth matrices. 
}

\references{
  Bowman, A. & Azzalini, A. (1997) \emph{Applied Smoothing Techniques
    for Data Analysis}. Oxford University Press, Oxford.

  Chacon, J.E. & Duong, T. (2018) \emph{Multivariate Kernel Smoothing
  and Its Applications}. Chapman & Hall/CRC, Boca Raton. 
  
  Duong, T. (2004) \emph{Bandwidth Matrices for Multivariate Kernel Density 
    Estimation.} Ph.D. Thesis, University of Western Australia. 

  Scott, D.W. (2015) \emph{Multivariate Density Estimation: Theory,
    Practice, and Visualization (2nd edn)}. John Wiley & Sons, New York.

  Silverman, B. (1986) \emph{Density Estimation for Statistics and
  Data Analysis}. Chapman & Hall/CRC, London.

  Simonoff, J. S. (1996) \emph{Smoothing Methods in Statistics}.
  Springer-Verlag, New York.

  Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman &
  Hall/CRC, London.
}

\keyword{package}

\seealso{\pkg{feature}, \pkg{sm}, \pkg{KernSmooth}}
