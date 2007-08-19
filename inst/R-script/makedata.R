### this script generates the example data sets
### (only for internal use)
### 19AUG2007: updated to "justvsn" but did not run code again.

rm(list=ls())
library(twilight)

### Leukemia data set of Golub et al. (1999)
library(golubEsets)
data(Golub_Merge)

### Variance-stabilizing normalization of Huber et al. (2002)
library(vsn)
golubNorm <- justvsn(Golub_Merge)

### A vector of class labels.
id <- as.numeric(Golub_Merge$ALL.AML)

expval <- twilight.pval(golubNorm,id)
save(expval,file="../../data/expval.rda")

exfdr <- twilight(expval,B=1000)
save(exfdr,file="../../data/exfdr.rda")
