\name{twilight.pval}
\alias{twilight.pval}
\title{ A function to compute p-values from expression sets }
\description{
A function to compute two-sample t, Z and fold change equivalent test statistics (paired or unpaired) and correlation coefficients. Based on permutations, expected test statistics as given in Tusher et al. (2001) and empirical \eqn{p}-values are computed. Additional output are \eqn{q}-values computed as given in Storey and Tibshirani (2003). The resulting object is of class \code{twilight} and can be passed to functions \code{twilight} or \code{plot.twilight}.
}
\usage{
twilight.pval(xin, yin, method = "fc", paired = FALSE, B = 1000, yperm = NULL, balance = FALSE, quant.ci = 0.95, s0=NULL, verbose = TRUE, filtering = FALSE)
}
\arguments{
  \item{xin}{ Either an expression set (\code{exprSet}) or a data matrix with rows corresponding to features and columns corresponding to samples. }
  \item{yin}{ A numerical vector containing class labels. The higher label denotes the case, the lower label the control samples to test case vs. control. For correlation scores, \code{yin} can be any numerical vector of length equal to the number of samples. }
  \item{method}{ Character string: \code{"fc"} for fold change equivalent test (that is log ratio test), \code{"t"} for t-test, and \code{"z"} for Z-test. With \code{"pearson"} or \code{"spearman"}, the test statistic is either Pearson's correlation coefficient or Spearman's rank correlation coefficient.}
  \item{paired}{ Logical value. Depends on whether the samples are paired. Ignored if \code{method="pearson"} or \code{method="spearman"}.}
  \item{B}{ Numerical value specifying the number of permutations. }
  \item{yperm}{ Optional matrix containing in each row a permutation of the class labels in binary(!) format for two-sample testing. For computation of correlation scores, the rows of \code{yperm} have to contain the appropriate values or ranks. Use this argument carefully! If \code{yperm} is specified, no other permutation will be done. Please note that the first row of \code{yperm} MUST be the input vector \code{yin}. Otherwise, the \eqn{p}-value calculation will be incorrect.}
  \item{balance}{ Logical value. Depends on whether balanced or unbalanced permutations should be done. Ignored if \code{method="pearson"} or \code{method="spearman"}.}
  \item{quant.ci}{ Probability value for confidence lines. Lines are symmetric and denote the \code{quant.ci}-quantile of maximal absolute differences between each permutatin and the expected scores. }
  \item{s0}{ Fudge factor for variance correction in the Z-test. Takes effect only if \code{method="z"}. If \code{s0=NULL}: The fudge factor is set to the median of the pooled standard deviations.}
  \item{verbose}{ Logical value for message printing. }
  \item{filtering}{ Logical value for filtering for permutations of class labels that produce a set of complete null scores. Invokes function \code{twilight.filtering}. If \code{yperm} is specified, no filtering will be done. Note that the filtering is done on unbalanced permutations even if \code{balance=TRUE}. }
}
\details{
See vignette.
}
\value{
Returns a \code{twilight} object consisting of a \code{data.frame} named \code{result} with variables
  \item{observed}{Observed test statistics.}
  \item{expected}{Mean of order statistics of the permutation statistics.}
  \item{candidate}{Binary vector. \code{"1"} for observations exceeding the confidence lines.}
  \item{pvalue}{Empirical \eqn{p}-values from two-sided hypothesis tests.}
  \item{qvalue}{\eqn{q}-values computed as described in Storey and Tibshirani (2003).}
  \item{index}{Index of the original ordering.}
Values are sorted by absolute \code{observed} scores.
}
\note{
Additional output consists of

\tabular{ll}{
  \code{ci.line} \tab Quantile corresponding to \code{quant.ci}, passed for plotting. \cr
  \code{pi0} \tab Estimated prior probability. \cr
  \code{call} \tab Character string of function arguments. \cr
  \code{quant.ci} \tab Passes \code{quant.ci} for plotting. \cr
}

The remaining slots are left free for function \code{twilight}.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Storey JD and Tibshirani R (2003): Statistical significance for genomewide studies, \emph{PNAS} \bold{100(16)}, 9440--9445.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, \emph{PNAS} \bold{98(9)}, 5116--5121.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight}}, \code{\link{plot.twilight}}, \code{\link{twilight.combi}}, \code{\link{twilight.filtering}}, \code{\link{expval}} }

\examples{
data(expval)  # twilight object created by twilight.pval on Golub data set
expval
plot(expval,which="scores")
plot(expval,which="qvalues")
}
\keyword{ nonparametric }
