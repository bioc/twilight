\name{twilight.getmatrix}
\alias{twilight.getmatrix}
\title{ Internal function to extract the gene expression matrix from an expression set }
\description{If the input object is of class \code{ExpressionSet}, the function extracts the gene expression matrix.}
\usage{
twilight.getmatrix(xin)
}
\arguments{
  \item{xin}{Input object of class \code{ExpressionSet}.}
}
\value{Returns the matrix of gene expression values.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}} }
\keyword{ internal }
