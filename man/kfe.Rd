\name{kfe}
\alias{kfe}
\alias{Hpi.kfe}
\alias{Hpi.diag.kfe}
\alias{hpi.kfe}
\title{Kernel functional estimate}

\description{
  Kernel functional estimate for 1- to 6-dimensional data.}

\usage{
kfe(x, G, deriv.order, inc=1, binned, bin.par, bgridsize, deriv.vec=TRUE,
    add.index=TRUE, verbose=FALSE)
Hpi.kfe(x, nstage=2, pilot, pre="sphere", Hstart, binned=FALSE, 
    bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
Hpi.diag.kfe(x, nstage=2, pilot, pre="scale", Hstart, binned=FALSE,
    bgridsize, amise=FALSE, deriv.order=0, verbose=FALSE, optim.fun="optim")
hpi.kfe(x, nstage=2, binned=FALSE, bgridsize, amise=FALSE, deriv.order=0)
}

\arguments{
  \item{x}{vector/matrix of data values}
  \item{nstage}{number of stages in the plug-in bandwidth selector (1 or 2)}
  \item{pilot}{"dscalar" = single pilot bandwidth (default) \cr
    "dunconstr" = single unconstrained pilot bandwidth}
  \item{pre}{"scale" = \code{\link{pre.scale}}, "sphere" = \code{\link{pre.sphere}}}
  \item{Hstart}{initial bandwidth matrix, used in numerical optimisation}
  \item{binned}{flag for binned estimation}
  \item{bgridsize}{vector of binning grid sizes}
  \item{amise}{flag to return the minimal scaled PI value}
  \item{deriv.order}{derivative order}
  \item{verbose}{flag to print out progress information. Default is FALSE.}
  \item{optim.fun}{optimiser function: one of \code{nlm} or \code{optim}}
  \item{G}{pilot bandwidth matrix}
  \item{inc}{0=exclude diagonal, 1=include diagonal terms in kfe
    calculation}
  \item{bin.par}{binning parameters - output from \code{\link{binning}}}
  \item{deriv.vec}{flag to compute duplicated partial derivatives in the 
    vectorised form. Default is FALSE.}
  \item{add.index}{flag to output derivative indices matrix. Default is true.}  
}

\value{
   Plug-in bandwidth matrix for \eqn{r}{r}-th order kernel functional estimator.
}

\details{
  \code{Hpi.kfe} is the optimal plug-in bandwidth for \eqn{r}{r}-th order kernel functional estimator
  based on the unconstrained pilot selectors of Chacon & Duong (2010).
  \code{hpi.kfe} is the 1-d equivalent, using the formulas from
  Wand & Jones (1995, p.70).

  \code{kfe} does not usually need to be called explicitly by the user.     
}

\references{
  Chacon, J.E. & Duong, T. (2010) Multivariate plug-in bandwidth
  selection with unconstrained pilot matrices. \emph{Test}, \bold{19},
  375-398.
  
  Wand, M.P. & Jones, M.C. (1995) \emph{Kernel Smoothing}. Chapman &
  Hall/CRC, London. 
}

\seealso{\code{\link{kde.test}}}
 
\keyword{smooth}
