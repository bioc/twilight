\name{twilight.permute.pair}
\alias{twilight.permute.pair}
\title{A permutation matrix of paired class labels}
\description{The function returns a matrix where each row is a (un)balanced permutation of the input twosample class labels.
}
\usage{
twilight.permute.pair(v, m, bal = TRUE)
}
\arguments{
  \item{v}{ A binary vector representing class labels in original order. Pairs must be in the same order. }
  \item{m}{ A numerical value giving the number of permutations. }
  \item{bal}{ Logical value. Results in balanced or unbalanced permutations. }
}
\value{Returns a matrix where each row contains one permuted vector of class labels. Note that even for balanced permutations, the first row always contains the original vector.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.permute.unpair}}, \code{\link{twilight.combi}} }
\keyword{ datagen }
