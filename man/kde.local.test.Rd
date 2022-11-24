\name{kde.local.test}
\alias{kde.local.test}
\title{Kernel density based local two-sample comparison test}

\description{
  Kernel density based local two-sample comparison test for 1- to 6-dimensional data.}

\usage{
kde.local.test(x1, x2, H1, H2, h1, h2, fhat1, fhat2, gridsize, binned, 
   bgridsize, verbose=FALSE, supp=3.7, mean.adj=FALSE, signif.level=0.05,
   min.ESS, xmin, xmax)
}

\arguments{
  \item{x1,x2}{vector/matrix of data values}
  \item{H1,H2,h1,h2}{bandwidth matrices/scalar bandwidths.  If these are missing, \code{Hpi} or \code{hpi} is called by default.}
  \item{fhat1,fhat2}{objects of class \code{kde}}
  \item{binned}{flag for binned estimation}
  \item{gridsize}{vector of grid sizes}
  \item{bgridsize}{vector of binning grid sizes}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{supp}{effective support for normal kernel}
  \item{mean.adj}{flag to compute second order correction for mean value of critical sampling distribution. Default is FALSE. Currently implemented for d<=2 only.}
  \item{signif.level}{significance level. Default is 0.05.}
  \item{min.ESS}{minimum effective sample size. See below for details.}
  \item{xmin,xmax}{vector of minimum/maximum values for grid} 
}

\value{
  A kernel two-sample local significance is an object of class
  \code{kde.loctest} which is a list with fields:
  \item{fhat1,fhat2}{kernel density estimates, objects of class \code{kde}}
  \item{chisq}{chi squared test statistic}
  \item{pvalue}{matrix of local \eqn{p}{p}-values at each grid point}
  \item{fhat.diff}{difference of KDEs}
  \item{mean.fhat.diff}{mean of the test statistic}  
  \item{var.fhat.diff}{variance of the test statistic}
  \item{fhat.diff.pos}{binary matrix to indicate locally significant fhat1 > fhat2}
  \item{fhat.diff.neg}{binary matrix to indicate locally significant fhat1 < fhat2}
  \item{n1,n2}{sample sizes}
  \item{H1,H2,h1,h2}{bandwidth matrices/scalar bandwidths}
}

\details{The null hypothesis is \eqn{H_0(\bold{x}): f_1(\bold{x}) = f_2(\bold{x})}{H_0(x): f_1(x) = f_2(x)} where \eqn{f_1, f_2}{f_1, f_2} 
  are the respective density functions. The measure of discrepancy is  
  \eqn{U(\bold{x}) = [f_1(\bold{x}) - f_2(\bold{x})]^2}{U(x) = [f_1(x) - f_2(x)]^2}. 
  Duong (2013) shows that the test statistic obtained, by substituting the
  KDEs for the true densities, has a null 
  distribution which is asymptotically chi-squared with 1 d.f.

  The required input is either \code{x1,x2} and \code{H1,H2}, or
  \code{fhat1,fhat2}, i.e. the data values and bandwidths or objects of class
  \code{kde}. In the former case, the \code{kde} objects are created.
  If the \code{H1,H2} are missing then the default are the plug-in
  selectors \code{Hpi}.  Likewise for missing \code{h1,h2}. 

  The \code{mean.adj} flag determines whether the
second order correction to the mean value of the test statistic should be computed. 
  \code{min.ESS} is borrowed from Godtliebsen et al. (2002)
  to reduce spurious significant results in the tails, though by it is usually
  not required for small to moderate sample sizes. 
}

\references{
  Duong, T. (2013) Local significant differences from non-parametric
  two-sample tests. \emph{Journal of Nonparametric Statistics},
  \bold{25}, 635-645.
    
  Godtliebsen, F., Marron, J.S. & Chaudhuri, P. (2002) 
  Significance in scale space for bivariate density estimation.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{11}, 1-22.

}

\seealso{\code{\link{kde.test}}, \code{\link{plot.kde.loctest}}}
  
\examples{
data(crabs, package="MASS")
x1 <- crabs[crabs$sp=="B", 4]
x2 <- crabs[crabs$sp=="O", 4]
loct <- kde.local.test(x1=x1, x2=x2)
plot(loct, ylim=c(-0.08,0.12))
cols <- hcl.colors(palette="Dark2",2)
plot(loct$fhat1, add=TRUE, col=cols[1])
plot(loct$fhat2, add=TRUE, col=cols[2])

## see examples in ? plot.kde.loctest
}

\keyword{test}
