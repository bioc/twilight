twilight.getlambda <- function(xin,verbose=TRUE){
### Function finds suitable regularization parameter "lambda". 
### For a sequence of lambdas, the objective function of SEP is
### computed for subsamples of the p-value vector "xin". 
### The final estimate is chosen based upon a Wilcoxon test
### comparison between objective function values of lambda=0 and
### each lambda>0. The penalized objective function should not
### differ a lot from the unpenalized one. Therefore, the 
### highest lambda that leads to a non-significant difference
### in means is chosen. The Wilcoxon ranksum test p-values are
### Bonferroni-corrected.

  ### Calling SEP.
  funk <- function(a,b){
    x   <- .C("sep",
              as.double(a),
              as.integer(length(a)),
              as.double(b),
              e = integer(length(a)),
              f = double(1),PACKAGE="twilight")$f
    return(x)
  }

  ### Try some values.
  lam  <- c(seq(0.005,0.195,by=0.005),seq(0.2,0.5,by=0.01))
  n    <- min(100,length(xin))
  
  boot.p <- matrix(NA,n,500)
  for (i in 1:500){
    boot.p[,i] <- sort(sample(xin,n,replace=TRUE))
  }

  obj0 <- apply(boot.p,2,funk,0)

  if (verbose){
    cat("Wait for",length(lam),"dots. \n")
  }

  wilc <- numeric(length(lam))
  for (i in 1:length(lam)){      
    obj.fun <- apply(boot.p,2,funk,lam[i])    
    wilc[i] <- wilcox.test(x=obj0,y=obj.fun)$p.value
    if (verbose){cat(".")}
  }
  if (verbose){cat("\n")}

  
  ### Bonferroni correction, look for non-significant differences.
  y <- which((wilc*length(wilc))<=0.05)

  if (length(y)==0){
    res <- NaN
  }
  if (length(y)!=0){
    res <- lam[y[1]-1]
  }
  
  return(res)
}
