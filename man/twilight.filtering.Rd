\name{twilight.filtering}
\alias{twilight.filtering}
\title{ Permutation filtering }
\description{
The function call invokes the filtering for permutations of class labels that produce a set of complete null scores. Depending on the test setting, the algorithm iteratively generates valid permutations of the class labels and computes scores. These are transformed to pooled p-values and each set of permutation p-values is tested for uniformity. Permutations with acceptable uniform p-value distributions are kept. The search stops if either the number \code{num.perm} of wanted permutations is reached or if the number of possible unique(!) permutations is smaller than \code{num.perm}. The default values are similar to function \code{twilight.pval} but please note the details below.
}
\usage{
twilight.filtering(xin, yin, method = "fc", paired = FALSE, s0 = 0, verbose = TRUE, num.perm = 1000, num.take = 50)
}
\arguments{
  \item{xin}{ Either an expression set (\code{ExpressionSet}) or a data matrix with rows corresponding to features and columns corresponding to samples. }
  \item{yin}{ A numerical vector containing class labels. The higher label denotes the case, the lower label the control samples to test case vs. control. For correlation scores, \code{yin} can be any numerical vector of length equal to the number of samples. }
  \item{method}{ Character string: \code{"fc"} for fold change equivalent test (that is log ratio test), \code{"t"} for t-test, and \code{"z"} for Z-test. With \code{"pearson"} or \code{"spearman"}, the test statistic is either Pearson's correlation coefficient or Spearman's rank correlation coefficient.}
  \item{paired}{ Logical value. Depends on whether the samples are paired. Ignored if \code{method="pearson"} or \code{method="spearman"}.}
 \item{s0}{ Fudge factor for variance correction in the Z-test. Takes effect only if \code{method="z"}. If \code{s0=0}: The fudge factor is set to the median of the pooled standard deviations.}
  \item{verbose}{ Logical value for message printing. }
  \item{num.perm}{ Number of permutations. Within \code{twilight.pval}, \code{num.perm} is set to \code{B}. }
  \item{num.take}{ Number of permutations kept in each step of the iterative filtering. Within \code{twilight.pval}, \code{num.take} is set to the minimum of 50 and \code{ceiling(num.perm/20)}. }
}
\details{
See vignette.
}
\value{
 Returns a matrix with permuted input labels \code{yin} as rows. Please note that this matrix is already translated into binary labels for two-sample testing or to ranks if Spearman's correlation was chosen. The resulting permutation matrix can be directly passed into function \code{twilight.pval}. Please note that the first row always contains the original input \code{yin} to be consistent with the other permutation functions in package \code{twilight}.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}} }
\examples{
\dontrun{
### Leukemia data set of Golub et al. (1999).
library(golubEsets)
data(Golub_Train)

### Variance-stabilizing normalization of Huber et al. (2002).
library(vsn)
golubNorm <- justvsn(Golub_Train)

### A binary vector of class labels.
id <- as.numeric(Golub_Train$ALL.AML)

### Do an unpaired t-test.
### Let's have a quick example with 50 filtered permutations only.
### With num.take=10, we only need 5 iteration steps.
yperm <- twilight.filtering(golubNorm,id,method="t",num.perm=50,num.take=10)
dim(yperm)

### Let's see if the filtered permutations really produce uniform p-value distributions.
### The first row is the original labeling, so we try the second permutation.
yperm <- yperm[-1,]
b <- twilight.pval(golubNorm,yperm[1,],method="t",yperm=yperm)
hist(b$result$pvalue)
}
}
\keyword{ datagen }
