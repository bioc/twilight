\name{expval}
\alias{expval}
\docType{data}
\title{ Example of twilight.pval result. }
\description{
Application of function twilight.pval on leukemia data set of Golub et al. (1999), as given in \code{data(golubMerge)} in \code{library(golubEsets)}.

First step was the variance-stabilizing normalization of Huber et al. (2002) in \code{library(vsn)}: \code{golubNorm <- vsn(exprs(golubMerge))}.

The function call was then \code{expval <- twilight.pval(golubNorm,id)} with \code{id <- as.numeric(golubMerge$ALL.AML)}.
}
\usage{data(expval)}
\format{A twilight object.}
\references{
Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, Coller H, Loh ML, Downing JR, Caligiuri MA, Bloomfield CD and Lander ES (1999): Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring, Science, vol. 286, pp. 531-537.

Huber W, von Heydebreck A, Sültmann H, Poustka A and Vingron M (2002): Variance stabilization applied to microarray data calibration and to the quantification of differential expression, Bioinformatics, vol. 18, suppl. 1, pp. S96-S104.
}
\keyword{datasets}
