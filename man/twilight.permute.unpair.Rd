\name{twilight.permute.unpair}
\alias{twilight.permute.unpair}
\title{A permutation matrix of unpaired class labels}
\description{The function returns a matrix where each row is a (un)balanced permutation of the input twosample class labels.
}
\usage{
twilight.permute.unpair(v, m, bal = TRUE)
}
\arguments{
  \item{v}{ A binary vector representing class labels in original order. }
  \item{m}{ A numerical value giving the number of permutations. }
  \item{bal}{ Logical value. Results in balanced or unbalanced permutations. }
}
\value{Returns a matrix where each row contains one permuted vector of class labels. Note that even for balanced permutations, the first row always contains the original vector.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.permute.pair}}, \code{\link{twilight.combi}} }
\keyword{ datagen }
