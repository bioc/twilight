\name{plot.twilight}
\alias{plot.twilight}
\title{A plot function for twilight objects.}
\description{Interface to plotting of \code{twilight} objects. Produces one of five possible plots.}
\usage{plot.twilight(x, which = "plot3", grayscale = FALSE, ...)}
\arguments{
  \item{x}{ An object of class twilight. }
  \item{which}{ A character string specifying the plot to be made. }
  \item{grayscale}{ Logical value. Specifying whether plots should be colored or grayscaled. Only necessary for "plot1" and "plot3". }
  \item{...}{ Additional graphical arguments. }
}
\details{
Option \code{which="plot1"} plots the expected versus the observed test statistics and draws confidence lines calculated from permutations. This plot is similar to plots in Tusher et al. (2001).

Option \code{which="plot2"} plots \eqn{q}-values versus the number of rejected hypotheses.

Option \code{which="plot3"} plots \eqn{u}-values versus 1 - local false discovery rate, that is the conditional probability of being significant given the corresponding \eqn{u}-value, plus bootstrap estimates if computed. Bottom ticks are 1% quantiles of \eqn{u}-values.

Option \code{which="plot4"} results in the volcano plot, that is observed score versus 1 - local false discovery rate. Bottom ticks are 1% quantiles of scores.

Option \code{which="plot5"} plots the observed fold change equivalent score distribution overlaid by the estimated effect size distribution, that is distribution of scores under the alternative. This plot is only available if function \code{twilight.pval} was applied with \code{method="fc"} as fold change equivalent scores are computed from log ratios.

Option \code{which="table"} tabulates histogram \code{"plot5"}.
}
\value{No value is returned.}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, submitted.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, PNAS 98(9), pp. 5116-5121.
 }
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}}, \code{\link{twilight}} }
\examples{
data(exfdr) # contains a twilight object created by twilight
plot(exfdr,which="plot3")
}
\keyword{ hplot }
