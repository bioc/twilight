\name{twilight.getmatrix}
\alias{twilight.getmatrix}
\title{Internal function to extract the gene expression matrix from an expression set}
\description{If the input object is of class \code{exprSet}, the function extracts the gene expression matrix.}
\usage{
twilight.getmatrix(xin)
}
\arguments{
  \item{xin}{Input object of class \code{exprSet}.}
}
\value{Returns the matrix of gene expression values.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}} }
\keyword{ internal }
