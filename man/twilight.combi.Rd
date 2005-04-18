\name{twilight.combi}
\alias{twilight.combi}
\title{ All permutations of a binary vector. }
\description{
For a given binary input vector, the function completely enumerates all possible permutations.
}
\usage{
twilight.combi(xin, pin, bin)
}
\arguments{
  \item{xin}{ Binary input vector, e.g. class labels. }
  \item{pin}{ Logical value. \code{TRUE} if samples are paired, \code{FALSE} if not. }
  \item{bin}{ Logical value. \code{TRUE} if permutations should be balanced, \code{FALSE} if not. }
}
\details{
Please note, that the resulting permutations are always as "balanced" as possible. The balancing is done for the smaller subsample. If its sample size is odd, say 5, \code{twilight.combi} computes all permutations with 2 or 3 samples unchanged. In the paired case, the output matrix contains only one half of all permutations. The second half is simply 1-\code{output} which leads to the same absolute test statistics in a paired test.
}
\value{
Returns a matrix where each row contains one permuted vector. Note that even for balanced permutations, the first row always contains the original vector. If the number of rows exceeds 10000, \code{NULL} is returned.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, IEEE TCBB 1(3), pp. 98-108.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.permute.pair}}, \code{\link{twilight.permute.unpair}} }

\examples{
x <- c(rep(0,4),rep(1,3))
y <- twilight.combi(x,pin=FALSE,bin=FALSE)
}
\keyword{ datagen }
