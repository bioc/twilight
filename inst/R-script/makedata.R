### this script generates the example data sets
### (only for internal use)

rm(list=ls())
library(twilight)

### Leukemia data set of Golub et al. (1999)
library(golubEsets)
data(golubMerge)

### Variance-stabilizing normalization of Huber et al. (2002)
library(vsn)
golubNorm <- vsn(exprs(golubMerge))

### A vector of class labels.
id <- as.numeric(golubMerge$ALL.AML)

expval <- twilight.pval(golubNorm,id)
save(expval,file="../../data/expval.rda")

exfdr <- twilight(expval,B=1000)
save(exfdr,file="../../data/exfdr.rda")
