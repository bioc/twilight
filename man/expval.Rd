\name{expval}
\alias{expval}
\title{ Example of twilight.pval result }
\docType{data}
\description{
Application of function twilight.pval on leukemia data set of Golub et al. (1999), as given in \code{data(Golub_Merge)} in \code{library(golubEsets)}.

First step was the variance-stabilizing normalization of Huber et al. (2002) in \code{library(vsn)}: \code{golubNorm <- vsn(exprs(Golub_Merge))}.

The function call was then \code{expval <- twilight.pval(golubNorm,id)} with \code{id <- as.numeric(Golub_Merge$ALL.AML)}.
}
\usage{data(expval)}
\format{A twilight object.}
\references{
Golub TR, Slonim DK, Tamayo P, Huard C, Gaasenbeek M, Mesirov JP, Coller H, Loh ML, Downing JR, Caligiuri MA, Bloomfield CD and Lander ES (1999): Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring, \emph{Science} \bold{286}, 531--537.

Huber W, von Heydebreck A, Sültmann H, Poustka A and Vingron M (2002): Variance stabilization applied to microarray data calibration and to the quantification of differential expression, \emph{Bioinformatics} \bold{18}, suppl. 1, S96--S104.

Scheid S and Spang R (2004): A stochastic downhill search algorithm for estimating the local false discovery rate, \emph{IEEE TCBB} \bold{1(3)}, 98--108.

Scheid S and Spang R (2005): twilight; a Bioconductor package for estimating the local false discovery rate, \emph{Bioinformatics} \bold{21(12)}, 2921--2922. 

Scheid S and Spang R (2006): Permutation filtering: A novel concept for significance analysis of large-scale genomic data, in: Apostolico A, Guerra C, Istrail S, Pevzner P, and Waterman M (Eds.): \emph{Research in Computational Molecular Biology: 10th Annual International Conference, Proceedings of RECOMB 2006, Venice, Italy, April 2-5, 2006}. Lecture Notes in Computer Science vol. 3909, Springer, Heidelberg, pp. 338-347.

Tusher VG, Tibshirani R and Chu G (2001): Significance analysis of mircroarrays applied to the ionizing response, \emph{PNAS} \bold{98(9)}, 5116--5121.
}
\keyword{datasets}
