\name{kde.test}
\alias{kde.test}
\alias{Hpi.kfe}
\alias{hpi.kfe}

\title{Kernel density based two-sample comparison test}

\description{
  Kernel density based two-sample comparison test for 1- to 6-dimensional data.}

\usage{
kde.test(x1, x2, H1, H2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, 
    binned=FALSE, bgridsize, verbose=FALSE, pilot="dscalar")
Hpi.kfe(x, nstage=2, pilot="dscalar", pre="sphere", Hstart, binned=FALSE, 
    bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="nlm")
hpi.kfe(x, nstage=2, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0)
}

\arguments{
  \item{x,x1,x2}{vector/matrix of data values}
  \item{H1,H2,h1,h2}{bandwidth matrices/scalar bandwidths. If these are missing, \code{Hpi.kfe} or \code{hpi.kfe} is called by default.}
  \item{psi1,psi2}{zero-th order kernel functional estimates}
  \item{var.fhat1,var.fhat2}{sample variance of KDE estimates evaluated at x1, x2}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{\code{"dscalar"} = single pilot bandwidth \cr
    \code{"dunconstr"} = single unconstrained pilot bandwidth}
  \item{pre}{\code{"scale"} = \code{\link{pre.scale}}, \code{"sphere"} = \code{\link{pre.sphere}}}
  \item{Hstart}{initial bandwidth matrix, used in numerical optimisation}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{deriv.order}{derivative order of kfe (kernel functional estimate). 
     Only deriv.order=0 is currently implemented.}
  \item{optim.fun}{optimiser function: one of \code{\link{nlm}} or \code{\link{optim}}.}
}

\value{
A list with fields
  \item{Tstat}{T statistic}
  \item{zstat}{z statistic - normalised version of Tstat}
  \item{pvalue}{p-value of the double sided test}
  \item{mean,var}{mean and variance of null distribution}
  \item{var.fhat1,var.fhat2}{sample variances of KDE values evaluated at data points}
  \item{n1,n2}{sample sizes}
  \item{H1,H2}{bandwidth matrices}
  \item{psi1,psi12,psi21,psi2}{kernel functional estimates}
}



\details{--The null hypothesis is \eqn{H_0: f_1 \equiv f_2}{f_1 = f_2} where \eqn{f_1, f_2}{f_1, f_2} 
  are the respective density functions. The measure of discrepancy is the integrated \eqn{L_2}{L2} error (ISE)
  \eqn{T = \int [f_1(\bold{x}) - f_2(\bold{x})]^2 \, d \bold{x}}{int [ f_1(x) - f_2(x)]^2 dx}. If 
  we rewrite this as \eqn{T = \psi_1 - \psi_{12} - \psi_{21} + \psi_2}{T = psi_1 - psi_12 - psi_21 + psi_2} 
  where \eqn{\psi_{uv} = \int f_u (\bold{x}) f_v (\bold{x})  \, d \bold{x}}{psi_uv = int f_u(x) f_v(x) dx},
  then we can use kernel functional estimators. Duong et al. (2012) show that this test statistic has a null 
  distribution which is asymptotically normal, so no bootstrap resampling is required to compute an approximate
  p-value.   
  
  As of ks 1.8.8, \code{kde.test(,binned=TRUE)} invokes binned estimation for
  the computation of the bandwidth selectors, and not the test statistic and 
  p-value. 

  --\code{Hpi.kfe} is the optimal plug-in bandwidth for \eqn{r}{r}-th order kernel functional estimator
  based on the unconstrained pilot selectors of Chacon & Duong (2010). This is automatically
  called by \code{kde.test} to estimate the \eqn{\psi}{psi} functionals with \eqn{r=0}{r=0}. \code{hpi.kfe} is the 1-d equivalent, using the formulas from
 Wand & Jones (1995, p.70).   
  
}

\references{
  Chacon, J.E. & Duong, T. (2010) Multivariate plug-in bandwidth
    selection with unconstrained pilot matrices. \emph{Test}, \bold{19}, 375-398.

  Duong, T., Goud, B & Schauer, K. (2012) Closed-form density-based framework for automatic detection of cellular morphology changes. \emph{PNAS}. In press.

  Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman & Hall/CRC, London.
}
  
\examples{
## univariate example
set.seed(8192)
samp <- 1000
x <- rnorm.mixt(n=samp, mus=0, sigmas=1, props=1)
x2 <- rnorm.mixt(n=samp, mus=0, sigmas=1, props=1)
y <- rnorm.mixt(n=samp, mus=0.25, sigmas=1, props=1)

kde.test(x1=x, x2=y)$pvalue   ## reject H0: f1=f2
kde.test(x1=x, x2=x2)$pvalue  ## accept H0: f1=f2


## bivariate example
mus1 <- rbind(c(1,-1), c(-1,1))
Sigmas1 <- rbind(invvech(c(4/9, 4/15, 4/9)), invvech(c(4/9, 4/15, 4/9)))
props1 <- c(1,1)/2
mus2 <- rbind(c(1,-1), c(-1,1))
Sigmas2 <- rbind(invvech(c(4/9, 14/45, 4/9)), 4/9*diag(2))
props2 <- c(1,1)/2

set.seed(8192)
samp <- 1000
x <- rmvnorm.mixt(n=samp, mus=mus1, Sigmas=Sigmas1, props=props1)
x2 <- rmvnorm.mixt(n=samp, mus=mus1, Sigmas=Sigmas1, props=props1)
y <- rmvnorm.mixt(n=samp, mus=mus2, Sigmas=Sigmas2, props=props2)
kde.test(x1=x, x2=y)$pvalue    ## reject H0: f1=f2
kde.test(x1=x, x2=x2)$pvalue   ## accept H0: f1=f2

## From ks 1.8.8, binned estimation is applied only to bandwidth
## matrices, not to the p-value calculation
kde.test(x1=x, x2=y, binned=TRUE)$pvalue 
}

\keyword{ test }
