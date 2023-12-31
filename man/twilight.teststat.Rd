\name{twilight.teststat}
\alias{twilight.teststat}
\title{ Interface to the test statistics provided within 'twilight' }
\description{
A function to compute two-sample t, Z and fold change equivalent test statistics (paired or unpaired) and correlation coefficients.
}
\usage{
twilight.teststat(xin, yin, method = "fc", paired = FALSE, s0 = NULL)
}
\arguments{
  \item{xin}{ Either an expression set (\code{ExpressionSet}) or a data matrix with rows corresponding to features and columns corresponding to samples. }
  \item{yin}{ A numerical vector containing class labels. The higher label denotes the case, the lower label the control samples to test case vs. control. For correlation scores, \code{yin} can be any numerical vector of length equal to the number of samples. }
  \item{method}{ Character string: \code{"fc"} for fold change equivalent test (that is log ratio test), \code{"t"} for t-test, and \code{"z"} for Z-test of Tusher et al. (2001). With \code{"pearson"} or \code{"spearman"}, the test statistic is either Pearson's correlation coefficient or Spearman's rank correlation coefficient.}
  \item{paired}{ Logical value. Depends on whether the samples are paired. Ignored if \code{method="pearson"} or \code{method="spearman"}.}
  \item{s0}{ Fudge factor for variance correction in the Z-test. Takes effect only if \code{method="z"}. If \code{s0=NULL}: The fudge factor is set to the median of the pooled standard deviations.}
}
\details{
Please see vignette for detailed information.
}
\value{
Returns a list with two components: a numerical vector of observed test statistics \code{observed}. Each entry corresponds to one row of the input data matrix. Also, the estimated fudge factor \code{s0} is returned. In any other case except \code{method="z"}, \code{s0} is zero.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, \emph{PNAS} \bold{98(9)}, 5116--5121.
}
\author{ Stefanie Scheid }
\seealso{ \code{\link{twilight.pval}} }

\examples{
### Z-test on random values
M <- matrix(rnorm(20000),nrow=1000)
id <- c(rep(1,10),rep(0,10))
stat <- twilight.teststat(M,id,method="z")

### Pearson correlation
id <- 1:20
stat <- twilight.teststat(M,id,method="pearson")
}
\keyword{ nonparametric }
