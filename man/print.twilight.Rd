\name{print.twilight}
\alias{print.twilight}
\title{A print function for twilight objects.}
\description{
Extract and print information about a twilight object.
}
\usage{
print.twilight(x, ...)
}
\arguments{
  \item{yin}{Input object of class \code{twilight}.}
  \item{...}{ Additional printing arguments. }
}
\value{
No value is returned.
}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, IEEE TCBB 1(3), pp. 98-108.
}
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{plot.twilight}} }
\examples{
data(exfdr) # contains a twilight object created by twilight
exfdr
}
\keyword{ internal }
