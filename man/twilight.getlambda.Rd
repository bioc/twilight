\name{twilight.getlambda}
\alias{twilight.getlambda}
\title{Internal function to compute regularization parameter lambda.}
\description{The SEP method optionally includes a penalty term in order to estimate the percentage of non-significant differences/genes.}
\usage{twilight.getlambda(xin, verbose = TRUE)}
\arguments{
  \item{xin}{A numerical vector containing \eqn{p}-values.}
  \item{verbose}{Logical value for message printing.}
}
\details{The function finds a suitable regularization parameter \eqn{\lambda}. For a sequence of lambdas, the objective function of SEP is computed for subsamples of the \eqn{p-value} vector \code{xin}. The final estimate is chosen based upon a Wilcoxon test comparison between objective function values of \eqn{\lambda=0} and each \eqn{\lambda>0}. The penalized objective function should not differ a lot from the unpenalized one. Therefore, the highest lambda that leads to a non-significant difference in means is chosen. The Wilcoxon ranksum test \eqn{p}-values are Bonferroni-corrected.}
\value{Returns numerical value.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, submitted.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight}} }
\keyword{ internal }
