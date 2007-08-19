\name{twilight}
\alias{twilight}
\title{ Estimation of the local false discovery rate }
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
In short, the successive exclusion procedure divides the set of p-values into two parts. The first part is chosen such that it resembles a uniform distribution while containing as many p-values as possible. This set of p-values represents the set of p-values derived from non-induced genes. The height of the uninform distribution is a natural estimate for the mixture parameter pi0. The p-value not contained in the uniform part serve as representatives of p-values derived from induced genes. Their distribution is the basis of the final estimate of the local false discovery rate.

NOTE: Library(snow) has to be loaded manually. It is not loaded as part of 'suggests' or 'depends' because it is only available under UNIX. If twilight does not work with the current version of snow, please send a report.
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

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.

Storey JD and Tibshirani R (2003): Statistical significance for genomewide studies, \emph{PNAS} \bold{100(16)}, 9440--9445.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}}, \code{\link{plot.twilight}}, \code{\link{exfdr}} }

\examples{
### twilight object created with B=1000 bootstrap samples
data(exfdr) 
print(exfdr)
plot(exfdr)
}
\keyword{ nonparametric }
