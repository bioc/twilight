\name{twilight}
\alias{twilight}
\title{ Estimation of local false discovery rate }
\description{
The function performs the successive exclusion procedure (SEP) as described in Scheid and Spang (2004).
}
\usage{
twilight(xin, lambda = NULL, B = 0, boot.ci = 0.95, clus = NULL, verbose = TRUE)
}
\arguments{
  \item{xin}{ Numerical vector of \eqn{p}-values or a \code{twilight} object. }
  \item{lambda}{ Numerical value denoting the penalty factor. If not specified, the function searchs for an appropriate regularization parameter. }
  \item{B}{ Numerical value specifying the number of bootstrap samples. If not specified, no bootstrap estimates are calculated. }
  \item{boot.ci}{ Numerical value denoting the probability value for bootstrap confidence intervals of local false discovery rate and prior pi0. }
  \item{clus}{ A list or numerical value to be passed to \code{makeCluster(clus,...)} in \code{library(snow)}. If specified, bootstrapping is performed in parallel. No checks on \code{clus} are performed. Please make sure that \code{makeCluster(clus)} works properly in your environment. }
  \item{verbose}{ Logical value for message printing. }
}
\details{
See vignette.
}
\value{
Returns a \code{twilight} object consisting of a \code{data.frame} named \code{result} with variables
  \item{pvalue}{Sorted input vector.}
  \item{qvalue}{\eqn{q}-values computed as described in Storey and Tibshirani (2003) with new estimate \code{pi0}.}
  \item{index}{Index of the original ordering.}
  \item{fdr}{Local false discovery rate averaged over 10 runs of SEP.}
  \item{mean.fdr}{Bootstrap estimate of local false discovery rate.}
  \item{lower.fdr}{Lower \code{boot.ci}-bootstrap confidence bound.}
  \item{upper.fdr}{ Upper \code{boot.ci}-bootstrap confidence bound.}
Values are sorted by \code{pvalue}.
}
\note{
Additional output consists of

\tabular{ll}{
  \code{lambda} \tab Regularization parameter. \cr
  \code{pi0} \tab SEP estimate of prior probability. \cr
  \code{boot.pi0} \tab Bootstrap estimate and \code{boot.ci}-bootstrap confidence bounds. \cr
  \code{boot.ci} \tab Passes \code{boot.ci} for plotting. \cr
  \code{effect} \tab Histogram of effect size distributions averaged over 10 runs of SEP. \cr
}

If \code{xin} is of class \code{twilight}, the remaining slots are filled with corresponding input values. If \code{xin} is not of class \code{twilight}, these slots remain free.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Storey JD and Tibshirani R (2003): Statistical significance for genomewide studies, \emph{PNAS} \bold{100(16)}, 9440--9445.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}}, \code{\link{plot.twilight}}, \code{\link{exfdr}} }

\examples{
data(exfdr) # twilight object created by twilight with B=1000 bootstraps
exfdr
plot(exfdr,which="fdr")
}
\keyword{ nonparametric }
