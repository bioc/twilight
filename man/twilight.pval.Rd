\name{twilight.pval}
\alias{twilight.pval}
\title{ A function to compute p-values from expression sets. }
\description{
The functions implements twosample t, Z and fold change equivalent tests, paired or unpaired. Based on permutations, expected test statistics as given in Tusher et al. (2001) and empirical \eqn{p}-values are computed. Additional output are \eqn{q}-values computed as given in Storey and Tibshirani (2003). The resulting object is of class \code{twilight} and can be passed to functions \code{twilight} or \code{plot.twilight}.
}
\usage{
twilight.pval(xin, yin, method = "fc", paired = FALSE, B = 10000, yperm = NULL, balance = FALSE, quant.ci = 0.95, s0=NULL, verbose = TRUE)
}
\arguments{
  \item{xin}{ Either an expression set (\code{exprSet}) or a data matrix with rows corresponding to features and columns corresponding to samples. }
  \item{yin}{ A numerical vector containing class labels. The higher label denotes the case, the lower label the control samples to test case vs. control. }
  \item{method}{ Character string: \code{"fc"} for fold change equivalent test (that is log ratio test), \code{"t"} for t test, and \code{"z"} for Z test. }
  \item{paired}{ Logical value. Depends on whether the samples are paired. }
  \item{B}{ Numerical value specifying the number of permutations. }
  \item{yperm}{ Optional matrix containing in each row a permutation of the class labels in binary format. If \code{yperm} is specified, no other permutation will be done. }
  \item{balance}{ Logical value. Depends on whether balanced or unbalanced permutations should be done. }
  \item{quant.ci}{ Probability value for confidence lines. Lines are symmetric and denote the \code{quant.ci}-quantile of maximal absolute differences between each permutatin and the expected scores. }
  \item{s0}{ Fudge factor for variance correction in Z test. Takes effect only if \code{method="z"}. If \code{s0=NULL}: Fudge factor is set to median of root pooled variances.}
  \item{verbose}{ Logical value for message printing. }
}
\details{
See vignette.
}
\value{
Returns a \code{twilight} object consisting of a \code{data.frame} named \code{result} with variables
  \item{observed}{Observed test statistics.}
  \item{expected}{Mean of order statistics of the permutation statistics.}
  \item{candidate}{Binary vector. \code{"1"} for observations exceeding the confidence lines.}
  \item{pvalue}{Twosided test permutation \eqn{p}-values.}
  \item{qvalue}{\eqn{q}-values computed as described in Storey and Tibshirani (2003).}
Values are sorted by \code{pvalue}.

Additional output consists of
  \item{ci.line}{Quantile corresponding to \code{quant.ci}, passed for plotting.}
  \item{pi0}{Estimated prior probability.}
  \item{call}{Character string of function arguments.}
  \item{quant.ci}{Passes \code{quant.ci} for plotting.}

The remaining slots are left free for function \code{twilight}.
}
\references{
Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, Coller H, Loh ML, Downing JR, Caligiuri MA, Bloomfield CD and Lander ES (1999): Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring, Science, vol. 286, pp. 531-537.

Huber W, von Heydebreck A, Sültmann H, Poustka A and Vingron M (2002): Variance stabilization applied to microarray data calibration and to the quantification of differential expression, Bioinformatics, vol. 18, suppl. 1, pp. S96-S104.

Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, IEEE TCBB 1(3), pp. 98-108.

Storey JD and Tibshirani R (2003): Statistical significance for genomewide studies, PNAS 100(16), pp. 9440-9445.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, PNAS 98(9), pp. 5116-5121.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight}}, \code{\link{plot.twilight}}, \code{\link{twilight.combi}} }

\examples{
### Leukemia data set of Golub et al. (1999)
library(golubEsets)
data(golubMerge)

### Variance-stabilizing normalization of Huber et al. (2002)
library(vsn)
golubNorm <- vsn(exprs(golubMerge))

### A vector of class labels.
id <- as.numeric(golubMerge$ALL.AML)

### Only 1000 permutations to reduce computation time.
a <- twilight.pval(golubNorm,id,B=1000)
plot(a,which="scores")
plot(a,which="qvalues")
}
\keyword{ nonparametric }
