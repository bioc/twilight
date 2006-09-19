twilight <- function(xin,lambda=NULL,B=0,boot.ci=0.95,clus=NULL,verbose=TRUE){
### Function computes local false discovery rates for vector of p-values.
###
### INPUT
###
### "xin":      Vector of p-values or object of class "twilight".
### "lambda":   Regularization paramater. If not specified, function
###             "twilight.getlambda" is called.
### "B":        Number of bootstrap runs. Non-positive values of "B"
###             correspond to no bootstrapping.
### "boot.ci":  Probability value for bootstrap confidence intervals of
###             local FDR and prior pi0.
### "clus":     Use your cluster to compute bootstraps in parallel.
### "verbose":  TRUE or FALSE.
###
### OUTPUT
###
### Object of class "twilight" containing a data.frame with
### "pvalue":     Sorted input vector. 
### "qvalue":     q-values computed as described in Remark B of:
###               Storey JD and Tibshirani R (2003): Statistical significance
###               for genomewide studies, PNAS 100(16), pp. 9440-9445.
### "fdr":        Local false discovery rate averaged over 10 runs of SEP.
### "mean.fdr":   Bootstrap estimate of local false discovery rate.
### "lower.fdr":  Lower "boot.ci"-bootstrap confidence bound.
### "upper.fdr":  Upper "boot.ci"-bootstrap confidence bound.
###
### Additional output:
### "lambda":     Regularization parameter.
### "pi0":        Estimated prior probability.
### "boot.pi0":   Bootstrap estimate and "boot.ci"-bootstrap confidence
###               bounds.
### "boot.ci":    Passes "boot.ci".
### "effect":     Histogram of effect size distributions averaged over 10 runs of SEP.
###
### If "xin" is of class "twilight", the remaining slots are filled
### with corresponding input values. If "xin" is not of class 
### "twilight", these slots remain free.

  ### Preliminary checks.
  if (class(xin)=="twilight"){
    pval  <- xin$result$pvalue
    score <- NULL
    check <- unlist(strsplit(xin$call," "))
    if (check[2]=="fc."){
      score <- xin$result$observed
    }
  }
  if (class(xin)!="twilight"){
    pval  <- xin
    rows  <- rownames(xin)
    score <- NULL
    index <- 1:length(pval)
    if (is.null(rownames(xin))){
      rows <- 1:length(pval)
    }
    if (is.unsorted(pval)){
      ix   <- order(pval)
      pval <- pval[ix]
      rows <- rows[ix]
      index <- index[ix]
    }
  }

  if (is.vector(pval)==FALSE){
    stop("Input must be vector of p-values or of class twilight.")
  }

  if ((max(pval)>1)|(min(pval)<0)){
    stop("Input vector must contain p-values in [0,1].")
  }

  
  ### Find suitable regularization parameter.
  if (is.null(lambda)){
    if (verbose){
      cat("Find suitable regularization parameter. ")
    }
    lambda <- twilight.getlambda(pval,verbose)
  }

  ### Compute the final estimate for pi0.
  funk.wrap <- function(a,b,c,d){

    ### Calling SEP. Returns binary vector (1=included,0=excluded).
    funk.sep <- function(ax,bx){
      x   <- .C("sep",
                as.double(ax),
                as.integer(length(ax)),
                as.double(bx),
                e = integer(length(ax)),
                f = double(1),PACKAGE="twilight")$e
      return(x)
    }

    bin.vec <- funk.sep(a,b)
    sep.pi0 <- sum(bin.vec)/length(bin.vec)    
    res     <- list(sep.pi0=sep.pi0)
    
    if (is.null(d)==FALSE){
      br         <- hist(d,plot=FALSE,br=100)$breaks
      sep.effect <- hist(d[as.logical(1-bin.vec)],plot=FALSE,br=br)   
      res        <- list(sep.pi0=sep.pi0,sep.effect=sep.effect)
    }

    if (verbose){cat(".")}    
    return(res)
  }
  
  ### 10 runs of SEP on the full input vector.
  if (verbose){cat("Run SEP. ")}

  orig.pval <- matrix(NA,length(pval),10)
  for (i in 1:10){
    orig.pval[,i] <- pval
  }

  if (is.null(clus)==TRUE){      
    if (verbose){cat("Wait for 10 dots. \n")}

    sep.H0 <- apply(orig.pval,2,funk.wrap,lambda,pval,score)

    if (verbose){cat("\n")}
  }
  
  if (is.null(clus)==FALSE){
    if (verbose){cat("\n")}

    ### library(snow)
    cl <- makeCluster(clus)
    clusterEvalQ(cl,library(twilight))

    .twilight.lambda.xxx <<- lambda
    .twilight.pval.xxx   <<- pval
    .twilight.score.xxx  <<- score
    
    clusterExport(cl,c(".twilight.lambda.xxx",".twilight.pval.xxx",".twilight.score.xxx"))
    
    sep.H0 <- parCapply(cl,orig.pval,funk.wrap,.twilight.lambda.xxx,.twilight.pval.xxx,.twilight.score.xxx)
    
    rm(.twilight.score.xxx,pos=1)
  }
  
  sep.pi0    <- numeric(10)
  sep.effect <- NaN
  for (i in 1:10){
    sep.pi0[i] <- sep.H0[[i]]$sep.pi0

    if (is.null(score)==FALSE){
      if (i==1){sep.effect <- sep.H0[[i]]$sep.effect}
      if (i!=1){sep.effect$counts <- sep.effect$counts + sep.H0[[i]]$sep.effect$counts}
    }

  }

  br.br <- c(0.01,0.02,0.04,0.05,0.1)    
  for (i in 1:length(br.br)){

    br <- quantile(pval,seq(0,1,by=br.br[i]))
      
    histmix <- hist(pval,plot=FALSE,br=br)
      
    x <- histmix$mids
    y <- 1/histmix$density

    if (sum(y==0)==0 & sum(y==Inf)==0 & sum(is.nan(y))==0){break}
  }
      
  smooth <- try(smooth.spline(x,y,df=7,w=1/x),silent=TRUE)
  
  if (class(smooth)=="try-error"){
    stop("Twilight cannot run properly.\n The number of features (genes, transcripts etc.) might be too small.\n Also, problems occur if the number of unique p-values is too small (e.g. in case of small sample sizes). Consider computing theoretical p-values using t.test or similar functions.")
  }

  sep.pi0 <- mean(sep.pi0)
  sep.H0  <- sep.pi0*predict(smooth,pval)$y
  
  sep.H0[which(sep.H0<0)] <- 0
  sep.H0[which(sep.H0>1)] <- 1

  if (is.null(score)==FALSE){
    ### mean histogram counts
    sep.effect$counts <- sep.effect$counts/10

    ### remove unchanged elements
    sep.effect$intensities <- NULL
    sep.effect$density     <- NULL
    sep.effect$xname       <- NULL
  }
  
  ### Bootstrap estimates.
  if (B>0){
    if (verbose){cat("Run SEP on bootstrap samples. ")}
    
    funk.wrap <- function(a,b,c){
      
      ### Calling SEP. Returns binary vector (1=included,0=excluded).
      funk.sep <- function(ax,bx){
        x   <- .C("sep",
                  as.double(ax),
                  as.integer(length(ax)),
                  as.double(bx),
                  e = integer(length(ax)),
                  f = double(1),PACKAGE="twilight")$e
        return(x)
      }
      
      bin.vec <- funk.sep(a,b)
      
      br.br <- c(0.01,0.02,0.04,0.05,0.1)
      
      for (i in 1:length(br.br)){
        
        br <- quantile(a,seq(0,1,by=br.br[i]))
        
        histmix <- hist(a,plot=FALSE,br=br)
        
        x <- histmix$mids
        y <- 1/histmix$density
        
        if (sum(y==0)==0 & sum(y==Inf)==0 & sum(is.nan(y))==0){break}
      }

      smooth <- try(smooth.spline(x,y,df=7,w=1/x),silent=TRUE)

      if (class(smooth)=="try-error"){
        stop("Twilight cannot run properly.\n The number of features (genes, transcripts etc.) might be too small.\n Also, problems occur if the number of unique p-values is too small (e.g. in case of small sample sizes). Consider computing theoretical p-values using t.test or similar functions.")
      }
      
      sep.pi0 <- sum(bin.vec)/length(bin.vec)
      sep.H0  <- sep.pi0*predict(smooth,c)$y
      
      sep.H0[which(sep.H0<0)] <- 0
      sep.H0[which(sep.H0>1)] <- 1

      if (verbose){cat(".")}

      res <- list(sep.pi0=sep.pi0,sep.H0=sep.H0)
      return(res)
    }

    
    boot.pval <- matrix(NA,length(pval),B)
    for (i in 1:B){
      boot.pval[,i] <- sort(sample(pval,replace=TRUE))
    }
    
    if (is.null(clus)==TRUE){      
      if (verbose){cat("Wait for",B,"dots. \n")}
      b.H0 <- apply(boot.pval,2,funk.wrap,lambda,pval)
      if (verbose){cat("\n")}
    }
    
    if (is.null(clus)==FALSE){
      if (verbose){cat("\n")}
      b.H0 <- parCapply(cl,boot.pval,funk.wrap,.twilight.lambda.xxx,.twilight.pval.xxx)
    }

    b.pi0 <- numeric(B)
    for (i in 1:B){
      b.pi0[i] <- b.H0[[i]]$sep.pi0
      b.H0[[i]]$sep.pi0 <- NULL
    }

    b.H0 <- as.data.frame(b.H0)
    
    mean.pi0 <- mean(b.pi0)
    ci.pi0   <- quantile(b.pi0,c((1-boot.ci)/2,boot.ci+(1-boot.ci)/2))

    mean.H0  <- apply(b.H0,1,mean)
    ci.H0    <- apply(b.H0,1,quantile,c((1-boot.ci)/2,boot.ci+(1-boot.ci)/2))

    lower.H0 <- ci.H0[1,]
    upper.H0 <- ci.H0[2,]
  }
  
  if (B<=0){
    m <- length(pval)

    mean.H0  <- rep(NaN,m)
    lower.H0 <- rep(NaN,m)
    upper.H0 <- rep(NaN,m)
    mean.pi0 <- NaN
    ci.pi0   <- rep(NaN,2)
  }

  if (is.null(clus)==FALSE){
    rm(.twilight.pval.xxx,pos=1)
    rm(.twilight.lambda.xxx,pos=1)    
    stopCluster(cl)
  }
  
  ### calculate q-values as described in Remark B of:
  ###
  ### Storey JD and Tibshirani R (2003): Statistical significance for
  ### genomewide studies, PNAS 100(16), pp. 9440-9445.
  ###
  ### BUT: set pi0 to SEP estimate.
  ###
  if (verbose){cat("Compute q-values. \n")}
  m       <- length(pval)
  qval    <- numeric(m)
  rp      <- rank(pval)
  qval[m] <- min(1,sep.pi0*pval[m]*m/rp[m])
  for (i in (m-1):1){
    qval[i] <- min(sep.pi0*pval[i]*m/rp[i],qval[i+1])
  }

  ### Prepare output.
  if (class(xin)=="twilight"){
    res <- xin

    res$result$qvalue     <- qval
    res$result$fdr        <- sep.H0
    res$result$mean.fdr   <- mean.H0
    res$result$lower.fdr  <- lower.H0
    res$result$upper.fdr  <- upper.H0
    res$lambda            <- lambda
    res$pi0               <- sep.pi0
    res$boot.pi0          <- data.frame(pi0=mean.pi0,lower.pi0=ci.pi0[1],upper.pi0=ci.pi0[2],row.names="")
    res$boot.ci           <- boot.ci
    res$effect            <- sep.effect
  }
  
  if (class(xin)!="twilight"){
    res <- list(result=data.frame(
                  observed=rep(NaN,m),
                  expected=rep(NaN,m),
                  candidate=rep(NaN,m),
                  pvalue=pval,
                  qvalue=qval,
                  fdr=sep.H0,
                  mean.fdr=mean.H0,
                  lower.fdr=lower.H0,
                  upper.fdr=upper.H0,
                  index=index,
                  row.names=rows),
                ci.line=NaN,
                quant.ci=NaN,
                lambda=lambda,
                pi0=sep.pi0,
                boot.pi0=data.frame(pi0=mean.pi0,lower.pi0=ci.pi0[1],upper.pi0=ci.pi0[2],row.names=""),
                boot.ci=boot.ci,
                effect=sep.effect,
                call=NaN)
    class(res) <- "twilight"
  }
  
  return(res)
}
