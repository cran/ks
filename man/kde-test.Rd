\name{kde.test}
\alias{kde.test}
\alias{Hpi.kfe}

\title{Kernel density based two-sample comparison test}

\description{
  Kernel density based two-sample comparison test for 2- to 6-dimensional data.}

\usage{
kde.test(x1, x2, H1, H2, psi1, psi2, fhat1, fhat2, var.fhat1, var.fhat2, 
    double.loop=FALSE, binned=FALSE, bgridsize, verbose=FALSE, pre.scale=FALSE)
Hpi.kfe(x, nstage=2, Hstart, deriv.order=0, binned=FALSE, bgridsize, 
    double.loop=FALSE, amise=FALSE, verbose=FALSE)
}

\arguments{
  \item{x,x1,x2}{matrices of data values}
  \item{H1,H2}{bandwidth matrices. If these are missing, \code{Hpi.kfe} is called by default.}
  \item{fhat1,fhat2}{gridded \code{kde} objects - output from \code{\link{kde}} (for 2- or 3-d data)}
  \item{psi1,psi2}{zero-th order kernel functional estimates}
  \item{var.fhat1,var.fhat2}{sample variance of KDE estimates evaluated at x1, x2}
  \item{double.loop}{flag to iterate through double loop for unbinned estimation. Default is FALSE.}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{pre.scale}{flag to pre-scale combined data set. Default is FALSE.}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{Hstart}{initial bandwidth matrix, used in numerical
    optimisation}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{deriv.order}{derivative order of kfe (kernel functional estimate). 
     Only deriv.order=0 is currently implemented.}
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
  are the respective density functions. The measure of discrepancy is the integrated L2 error (ISE)
  \eqn{T = \int [f_1(\bold{x}) - f_2(\bold{x})]^2 \, d \bold{x}}{int [ f_1(x) - f_2(x)]^2 dx}. If 
  we rewrite this as \eqn{T = \psi_1 - \psi_{12} - \psi_{21} + \psi_2}{T = psi_1 - psi_12 - psi_21 + psi_2} 
  where \eqn{\psi_{uv} = \int f_u (\bold{x}) f_v (\bold{x})  \, d \bold{x}}{psi_uv = int f_u(x) f_v(x) dx},
  then we can use kernel functional estimators. Duong (2011) shows that this test statistic has a null 
  distribution which is asymptotically normal, so no bootstrap resampling is required to compute an approximate
  p-value.   

  --\code{Hpi.kfe} is the optimal plug-in bandwidth for \eqn{r}{r}-th order kernel functional estimator
  based on the unconstrained pilot selectors of Chacon & Duong (2010). This is automatically
  called by \code{kde.test} to estimate the \eqn{\psi}{psi} functionals with \eqn{r=0}{r=0}. 
}

\references{
  Chacon, J.E. & Duong, T. (2010) Multivariate plug-in bandwidth
    selection with unconstrained pilot matrices. \emph{Test}, \bold{19}, 375-398.

  Duong, T. (2011) Multivariate kernel density based two-sample comparisons without 
    bootstrap resampling. Submitted.
}
  
\examples{
mus1 <- rbind(c(1,-1), c(-1,1))
Sigmas1 <- rbind(invvech(c(4/9, 4/15, 4/9)), invvech(c(4/9, 4/15, 4/9)))
props1 <- c(1,1)/2
mus2 <- rbind(c(1,-1), c(-1,1))
Sigmas2 <- rbind(invvech(c(4/9, 14/45, 4/9)), 4/9*diag(2))
props2 <- c(1,1)/2

set.seed(8192)
samp <- 1000
x <- rmvnorm.mixt(n=samp, mus=mus1, Sigmas=Sigmas1, props=props1)
y <- rmvnorm.mixt(n=samp, mus=mus2, Sigmas=Sigmas2, props=props2)
y2 <- rmvnorm.mixt(n=samp, mus=mus2, Sigmas=Sigmas2, props=props2)
kde.test(x1=x, x2=y)$pvalue    ## reject H0: f1=f2
kde.test(x1=y2, x2=y)$pvalue   ## accept H0: f1=f2
}

\keyword{ test }
