twilight.pval <- function(xin,yin,method="fc",paired=FALSE,B=10000,yperm=NULL,balance=FALSE,quant.ci=0.95,s0=NULL,verbose=TRUE){ 
### Function computes test statistics for paired or unpaired twosided
### t-test, Z-test or Fold change test.
###
### INPUT
###
### "xin":      Data matrix.
### "yin":      Vector of real valued group labels.
###
### Higher valued group labels are used as case samples and lower
### valued group label as control samples to test Case vs. Control.
###
### "method":   "t" for t statistic, "z" for Z statistic and "fc" for
###             Fold change equivalent (that is log ratio).
### "paired":   TRUE or FALSE. Depends on test setting.
### "B":        Number of permutations to estimate the null distribution
###             and the expected test statistics.
### "yperm":    Optional matrix of permuted class labels 0 and 1. Each row must
###             contain one permutation. No other permutations will be done.
### "balance":  TRUE or FALSE. Gives balanced or unbalanced permutations.
### "quant.ci": Probability value for confidence lines. Lines are symmetric
###             and denote the "quant.ci"-quantile of maximal absolute
###             differences between each permutation and the expected scores.
### "s0":       Fudge factor for Z test. If s0=NULL, s0 is set to median of
###             root pooled variances in test statistic.
### "verbose":  TRUE or FALSE.
###
### OUTPUT
###
### Object of class "twilight" containing a data.frame with
### "observed":  Observed test statistics.
### "expected":  Mean of order statistics of the permutation statistics,
###              computed as described in:
###              Tusher VG, Tibshirani R and Chu G (2001): Significance
###              analysis of mircroarrays applied to the ionizing response,
###              PNAS 98(9), pp. 5116-5121.
### "candidate": Binary vector. "1" for genes exceeding the confidence lines.
### "pvalue":    Twosided test permutation p-values.
### "qvalue":    q-values computed as described in Remark B of:
###              Storey JD and Tibshirani R (2003): Statistical significance
###              for genomewide studies, PNAS 100(16), pp. 9440-9445.
###
### Additional output:
### "ci.line":   Quantile corresponding to "quant.ci".
### "pi0":       Estimated prior probability.
### "call":      String of function arguments.
### "quant.ci":  Passes "quant.ci".
###
### The remaining slots are left free for function "twilight".
  
  ### extract data matrix if class(xin) is exprSet
  xin <- twilight.getmatrix(xin)

  ### check dimensions
  if (is.matrix(xin)==FALSE){
    stop("First input must be a matrix. \n")
  }

  if (is.vector(yin)==FALSE){
    stop("Second input must be a vector. \n")
  }

  if (length(yin)!=dim(xin)[2]){
    if (length(yin)!=dim(xin)[1]){
      stop("Dimensions of input matrix and length of index vector do not match. \n")
    }    
    xin <- t(xin)
  }

  if (length(yin)!=dim(xin)[2]){
    if (length(yin)!=dim(xin)[1]){
      stop("Dimensions of input matrix and length of index vector do not match. \n")
    }    
  }

  if (is.null(yperm)==FALSE){
    if (length(yin)!=dim(yperm)[2]){
      stop("Dimensions of permutation matrix and length of index vector do not match. \n")
    }
    
    y <- unique(yperm[1,])  
    if ((length(y)!=2)|(min(y)!=0)|(max(y)!=1)){
      stop("Labels in permutation matrix must be 0 and 1. \n")
    }
  }
  
  ### translate index vector yin into binary vector with 1 as case and 0 as control samples.
  y <- sort(unique(yin),na.last=TRUE)

  if (length(y)!=2){
    stop("Samples must belong to TWO classes. \n")
  }

  y1 <- which(yin==y[1])
  y2 <- which(yin==y[2])
  yin[y1] <- 0
  yin[y2] <- 1

  if (paired==TRUE){
    if (sum(yin)!=sum(1-yin)){
      stop("This is a PAIRED twosample test. \n")
    }
  }

  ### prepare matrix of permuted index labels.
  if (is.null(yperm)){
    yperm <- twilight.combi(yin,pin=paired,bin=balance)
    if ((is.null(yperm)==FALSE)&(verbose)){cat("Complete enumeration possible. \n")}
  }
  if ((is.null(yperm))&(verbose)){cat("No complete enumeration. Prepare permutation matrix. \n")}
  if ((is.null(yperm))&(paired==FALSE)){yperm <- twilight.permute.unpair(yin,B,balance)}
  if ((is.null(yperm))&(paired==TRUE)){yperm <- twilight.permute.pair(yin,B,balance)}

  B <- dim(yperm)[1]

  ### compute observed test statistics.
  if (paired==FALSE){
    funk <- function(a,b,c,d,s){
      .C("unpaired",
         as.integer(a),
         as.integer(sum(a)),
         as.integer(sum(1-a)),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         as.integer(c),
         as.integer(which(d==1)-1),
         as.integer(which(d==0)-1),
         as.double(s),
         e=double(nrow(b)),PACKAGE="twilight")$e
    }
  }
  if (paired==TRUE){
    funk <- function(a,b,c,d,s){
      .C("paired",
         as.integer(a),
         as.integer(sum(a)),
         as.integer(sum(1-a)),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         as.integer(c),
         as.integer(which(d==1)-1),
         as.integer(which(d==0)-1),
         as.double(s),
         e=double(nrow(b)),PACKAGE="twilight")$e
    }
  }
  
  if (verbose){cat("Compute vector of observed statistics. \n")}

  ### Z test with s0=0 is a t test.
  if (is.null(s0)==FALSE){
    if((s0==0)&(method=="z")){method <- "t"}
  }
  
  if (is.null(s0)){s0 <- 0}
  stat.obs <- switch(method,
                     t = funk(yin,xin,1,yin,s0),
                     z = funk(yin,xin,2,yin,s0),
                    fc = funk(yin,xin,3,yin,s0))


  ### compute twosided test p-values from permutations.  
  ### sort permutation scores and calculate expected scores as described in:
  ### 
  ### Tusher VG, Tibshirani R and Chu G (2001): Significance
  ### analysis of mircroarrays applied to the ionizing response,
  ### PNAS 98(9), pp. 5116-5121.
  ###
  if (paired==FALSE){
    funk <- function(a,b,c,d,orig,s){
      x <- .C("unpairedperm",
              as.integer(t(a)),
              as.integer(nrow(a)),
              as.integer(sum(a[1,])),
              as.integer(sum(1-a[1,])),
              as.double(t(b)),
              as.integer(nrow(b)),
              as.integer(ncol(b)),
              as.integer(c),
              as.double(d),
              as.integer(which(orig==1)-1),
              as.integer(which(orig==0)-1),
              as.double(s),
              e=double(nrow(b)),
              f=double(nrow(b)),PACKAGE="twilight"
              )
      res <- list(exp=x$e,pval=x$f)
      return(res)
    }
  }
  if (paired==TRUE){
    funk <- function(a,b,c,d,orig,s){
      x <- .C("pairedperm",
              as.integer(t(a)),
              as.integer(nrow(a)),
              as.integer(sum(a[1,])),
              as.integer(sum(1-a[1,])),
              as.double(t(b)),
              as.integer(nrow(b)),
              as.integer(ncol(b)),
              as.integer(c),
              as.double(d),
              as.integer(which(orig==1)-1),
              as.integer(which(orig==0)-1),
              as.double(s),
              e=double(nrow(b)),
              f=double(nrow(b)),PACKAGE="twilight"
              )
      res <- list(exp=x$e,pval=x$f)
      return(res)
    }
  }
  
  if (verbose){cat("Compute expected scores and p-values. \n")}
  stat.exp <- switch(method,
                     t = funk(yperm,xin,1,stat.obs,yin,s0),
                     z = funk(yperm,xin,2,stat.obs,yin,s0),
                     fc = funk(yperm,xin,3,stat.obs,yin,s0))


  pval     <- stat.exp$pval
  stat.exp <- stat.exp$exp
  stat.exp <- stat.exp[rank(stat.obs)]
  
  ### sort all values according to the p-values.
  ix <- order(pval)
  stat.obs <- stat.obs[ix]
  stat.exp <- stat.exp[ix]
  pval     <- pval[ix]
  rows     <- rownames(xin)[ix]

  if (is.null(rownames(xin))){
    rows <- ix
  }
  
  ### calculate q-values as described in Remark B of:
  ###
  ### Storey JD and Tibshirani R (2003): Statistical significance for
  ### genomewide studies, PNAS 100(16), pp. 9440-9445.
  ###
  if (verbose){cat("Compute q-values. \n")}
  l <- seq(0,0.95,by=0.01)
  m <- length(pval)
  p <- numeric()
  for (i in 1:length(l)){
    p <- c(p,sum(pval>l[i])/(m*(1-l[i])))
  }  

  pi.model <- lm(p ~ ns(l,df=3))
  pi0 <- predict(pi.model,data.frame(l=c(l,1)))
  pi0 <- min(pi0[length(pi0)],1)

  qval    <- numeric(m)
  rp      <- rank(pval)
  qval[m] <- min(1,pi0*pval[m]*m/rp[m])
  for (i in (m-1):1){
    qval[i] <- min(pi0*pval[i]*m/rp[i],qval[i+1])
  }

  ### compute permutation based confidence lines for plot1.
  if (verbose){cat("Compute values for confidence lines. \n")}
  if (paired==FALSE){
    funk <- function(a,b,c,d,orig){
      x <- .C("unpairedci",
              as.integer(t(a)),
              as.integer(nrow(a)),
              as.integer(sum(a[1,])),
              as.integer(sum(1-a[1,])),
              as.double(t(b)),
              as.integer(nrow(b)),
              as.integer(ncol(b)),
              as.integer(c),
              as.double(d),
              as.integer(which(orig==1)-1),
              as.integer(which(orig==0)-1),
              e=double(nrow(a)),PACKAGE="twilight"
              )$e
    }
  }
  if (paired==TRUE){
    funk <- function(a,b,c,d,orig){
      x <- .C("pairedci",
              as.integer(t(a)),
              as.integer(nrow(a)),
              as.integer(sum(a[1,])),
              as.integer(sum(1-a[1,])),
              as.double(t(b)),
              as.integer(nrow(b)),
              as.integer(ncol(b)),
              as.integer(c),
              as.double(d),
              as.integer(which(orig==1)-1),
              as.integer(which(orig==0)-1),
              e=double(nrow(a)),PACKAGE="twilight"
              )$e
    }
  }

  ### compute confidence bounds
  ci.sel  <- sample(1:B,min(1000,B))
  ci.line <- switch(method,
                    t = funk(yperm[ci.sel,],xin,1,stat.exp,yin),
                    z = funk(yperm[ci.sel,],xin,2,stat.exp,yin),
                    fc = funk(yperm[ci.sel,],xin,3,stat.exp,yin))

  ci.line <- quantile(ci.line,quant.ci)

  ### mark genes with differences exceeding the confidence lines.
  cand <- as.numeric( abs(stat.obs-stat.exp)>ci.line )
  
  res <- list(result=data.frame(
                observed=stat.obs,
                expected=stat.exp,
                candidate=cand,
                pvalue=pval,
                qvalue=qval,
                fdr=rep(NaN,m),
                mean.fdr=rep(NaN,m),
                lower.fdr=rep(NaN,m),
                upper.fdr=rep(NaN,m),
                row.names=rows),
              ci.line=ci.line,
              quant.ci=quant.ci,
              lambda=NaN,
              pi0=pi0,
              boot.pi0=NaN,
              boot.ci=NaN,
              effect=NaN,
              call=paste("Test: ",method,". Paired: ",paired,". Number of permutations: ",B,". Balanced: ",balance,".",sep=""))
  class(res) <- "twilight"
  
  return(res)
}
