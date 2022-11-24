\name{kde.test}
\alias{kde.test}
\title{Kernel density based global two-sample comparison test}

\description{
  Kernel density based global two-sample comparison test for 1- to 6-dimensional data.}

\usage{
kde.test(x1, x2, H1, H2, h1, h2, psi1, psi2, var.fhat1, var.fhat2, 
    binned=FALSE, bgridsize, verbose=FALSE)
}

\arguments{
  \item{x1,x2}{vector/matrix of data values}
  \item{H1,H2,h1,h2}{bandwidth matrices/scalar bandwidths. If these are
  missing, \code{Hpi.kfe}, \code{hpi.kfe} is called by default.}
  \item{psi1,psi2}{zero-th order kernel functional estimates}
  \item{var.fhat1,var.fhat2}{sample variance of KDE estimates evaluated at x1, x2}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
}

\value{
A kernel two-sample global significance test is a list with fields:
  \item{Tstat}{T statistic}
  \item{zstat}{z statistic - normalised version of Tstat}
  \item{pvalue}{\eqn{p}{p}-value of the double sided test}
  \item{mean,var}{mean and variance of null distribution}
  \item{var.fhat1,var.fhat2}{sample variances of KDE values evaluated at data points}
  \item{n1,n2}{sample sizes}
  \item{H1,H2}{bandwidth matrices}
  \item{psi1,psi12,psi21,psi2}{kernel functional estimates}
}

\details{The null hypothesis is \eqn{H_0: f_1 \equiv f_2}{H_0: f_1 = f_2} where \eqn{f_1, f_2}{f_1, f_2} 
  are the respective density functions. The measure of discrepancy is
  the integrated squared error (ISE)
  \eqn{T = \int [f_1(\bold{x}) - f_2(\bold{x})]^2 \, d \bold{x}}{int [ f_1(x) - f_2(x)]^2 dx}. If 
  we rewrite this as \eqn{T = \psi_{0,1} - \psi_{0,12} - \psi_{0,21} + \psi_{0,2}}{T = psi_0,1 - psi_0,12 - psi_0,21 + psi_0,2} 
  where \eqn{\psi_{0,uv} = \int f_u (\bold{x}) f_v (\bold{x})  \, d \bold{x}}{psi_0,uv = int f_u(x) f_v(x) dx},
  then we can use kernel functional estimators. This test statistic has a null 
  distribution which is asymptotically normal, so no bootstrap
  resampling is required to compute an approximate \eqn{p}{p}-value.     
  
  If \code{H1,H2} are missing then the  plug-in selector \code{\link{Hpi.kfe}}
  is automatically called by \code{kde.test} to estimate the
  functionals with \code{kfe(, deriv.order=0)}. Likewise for missing
  \code{h1,h2}. 

  For \pkg{ks} \eqn{\geq}{>=} 1.8.8, \code{kde.test(,binned=TRUE)} invokes binned
  estimation for the computation of the bandwidth selectors, and not the
  test statistic and \eqn{p}{p}-value.  
}

\references{
  Duong, T., Goud, B. & Schauer, K. (2012) Closed-form density-based framework for automatic detection of cellular morphology changes. \emph{PNAS}, \bold{109}, 8382-8387. 
}

\seealso{\code{\link{kde.local.test}}}
 
\examples{
set.seed(8192)
samp <- 1000
x <- rnorm.mixt(n=samp, mus=0, sigmas=1, props=1)
y <- rnorm.mixt(n=samp, mus=0, sigmas=1, props=1)
kde.test(x1=x, x2=y)$pvalue   ## accept H0: f1=f2

library(MASS)
data(crabs)
x1 <- crabs[crabs$sp=="B", c(4,6)]
x2 <- crabs[crabs$sp=="O", c(4,6)]
kde.test(x1=x1, x2=x2)$pvalue  ## reject H0: f1=f2
}

\keyword{test}
