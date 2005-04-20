\name{plot.twilight}
\alias{plot.twilight}
\title{A plot function for twilight objects}
\description{Interface to plotting of \code{twilight} objects. Produces one of five possible plots.}
\usage{plot.twilight(x, which = "fdr", grayscale = FALSE, legend = TRUE, ...)}
\arguments{
  \item{x}{ An object of class twilight. }
  \item{which}{ A character string specifying the plot to be made. }
  \item{grayscale}{ Logical value. Specifying whether plots should be colored or grayscaled. Only necessary for "scores" and "fdr". }
  \item{legend}{ Logical value. Produces a legend for "scores" and "effectsize". A legend for "fdr" is only available if bootstrap estimates exist. }
  \item{...}{ Additional graphical arguments. }
}
\details{
Option \code{which="scores"} plots the expected versus the observed test statistics and draws confidence lines calculated from permutations. This plot is similar to plots in Tusher et al. (2001).

Option \code{which="qvalues"} plots \eqn{q}-values versus the number of rejected hypotheses.

Option \code{which="fdr"} plots \eqn{u}-values versus 1 - local false discovery rate, that is the conditional probability of being significant given the corresponding \eqn{u}-value, plus bootstrap estimates if computed. Bottom ticks are 1% quantiles of \eqn{u}-values.

Option \code{which="volcano"} results in the volcano plot, that is observed score versus 1 - local false discovery rate. Bottom ticks are 1% quantiles of scores.

Option \code{which="effectsize"} plots the observed fold change equivalent score distribution overlaid by the estimated effect size distribution, that is distribution of scores under the alternative. This plot is only available if function \code{twilight.pval} was applied with \code{method="fc"} as fold change equivalent scores are computed from log ratios.

Option \code{which="table"} tabulates histogram \code{"effectsize"}.
}
\value{No value is returned except for "table".}
\references{
Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, \emph{PNAS} \bold{98(9)}, 5116--5121.
 }
\author{ Stefanie Scheid \url{http://www.molgen.mpg.de/~scheid} }
\seealso{ \code{\link{twilight.pval}}, \code{\link{twilight}} }
\examples{
data(exfdr) # contains a twilight object created by twilight
plot(exfdr,which="fdr")
}
\keyword{ hplot }
