\name{kde.local.test}
\alias{kde.local.test}

\title{Kernel density based local two-sample comparison test}

\description{
  Kernel density based local two-sample comparison test for 1- to 6-dimensional data.}

\usage{
kde.local.test(x1, x2, H1, H2, h1, h2, fhat1, fhat2, binned=FALSE, 
   bgridsize, verbose=FALSE, supp=3.7, mean.adj=FALSE, 
   signif.level=0.05, min.ESS)
}

\arguments{
  \item{x1,x2}{vector/matrix of data values}
  \item{H1,H2,h1,h2}{bandwidth matrices/scalar bandwidths}
  \item{fhat1,fhat2}{objects of class \code{kde}}
  \item{binned}{flag for binned estimation. Default is FALSE.}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{supp}{effective support for normal kernel [-supp, supp]}
  \item{mean.adj}{flag to compute second order correction for mean value of critical sampling distribution. Default is FALSE. Currently implemented for d<=2 only.}
  \item{signif.level}{significance level. Default in 0.05.}
  \item{min.ESS}{minimum effective sample size. See below for details.}  
}

\value{
An object of class \code{kde.loctest} which is a list with fields
  \item{fhat1,fhat2}{kernel density estimates, objects of class \code{kde}}
  \item{chisq}{chi squared test statistic}
  \item{pvalue}{matrix of local p-values at each grid point}
  \item{fhat.diff}{difference of KDEs}
  \item{mean.fhat.diff}{mean of the test statistic}  
  \item{var.fhat.fiff}{variance of the test statistic}
  \item{fhat.diff.pos}{binary matrix to indicate locally signficant fhat1 > fhat2}
  \item{fhat.diff.neg}{binary matrix to indicate locally signficant fhat1 < fhat2}
  \item{n1,n2}{sample sizes}
  \item{H1,H2,h1,h2}{bandwidth matrices/bandwidths}
}



\details{The null hypothesis is \eqn{H_0(\bold{x}): f_1(\bold{x}) = f_2(\bold{x})}{H_0(x): f_1(x) = f_2(x)} where \eqn{f_1, f_2}{f_1, f_2} 
  are the respective density functions. The measure of discrepancy is  
  \eqn{U(\bold{x}) = [f_1(\bold{x}) - f_2(\bold{x})]^2}{U(x) = [f_1(x) - f_2(x)]^2}. 
  Duong (2012) show that the test statistic obtained, by substituting the
KDEs for the true densities, has a null 
  distribution which is asymptotically chi-squared with d d.f, where d is
  the data dimension.    

  The required input is either \code{x1,x2} and \code{H1,H2}, or \code{fhat1,fhat2}. 
  That is, the data values and bandwidths or a KDE object of class \code{kde}. 
  In the former case, a \code{kde} object is created with plug-in bandwidth
  \code{Hpi()}.  

  The \code{mean.adj} flag determines whether the
second order correction to the mean value of the test statistic should be computed. 
  \code{min.ESS} is borrowed from Godtliebsen et al. (2002)
  to reduce spurious significant results in the tails, though by it is usually
  not required for small to moderate sample sizes. 
}

\references{
  Duong, T. (2012) Local signficant differences from non-parametric
  two-sample tests. Submitted.
    
  Godtliebsen, F., Marron, J.S. & Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.

}

\seealso{\code{\link{kde.test}}}

  
\examples{
## univariate example
mus1 <- 0
sigmas1 <- 1
props1 <-1
    
mus2 <- c(0, -1, 1)
sigmas2 <- c(1, 1/4, 1/4)
props2 <- c(1/2, 1/4, 1/4)

x <- rnorm.mixt(n=1000, mus=mus1, sigmas=sigmas1, props=props1)
y <- rnorm.mixt(n=1000, mus=mus2, sigmas=sigmas2, props=props2)
loct <- kde.local.test(x1=x, x2=y, binned=TRUE)
plot(loct, lcol=2)     

## bivariate example
mus1 <- rbind(c(1,-1), c(-1,1))
Sigmas1 <- rbind(invvech(c(4/9, 4/15, 4/9)), invvech(c(4/9, 4/15, 4/9)))
props1 <- c(1,1)/2
    
mus2 <- rbind(c(1,-1), c(-1,1))
Sigmas2 <- rbind(invvech(c(4/9, 14/45, 4/9)), 4/9*diag(2))
props2 <- c(1,1)/2

x <- rmvnorm.mixt(n=10000, mus=mus1, Sigmas=Sigmas1, props=props1)
y <- rmvnorm.mixt(n=10000, mus=mus2, Sigmas=Sigmas2, props=props2)
loct <- kde.local.test(x1=x, x2=y, binned=TRUE)
plot(loct)
}

\keyword{ test }
