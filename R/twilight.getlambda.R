twilight.getlambda <- function(xin,verbose=TRUE){
### Function finds suitable regularization parameter "lambda". 
### For a sequence of lambdas, the objective function of SEP is
### computed for subsamples of the p-value vector "xin". 
### The final estimate is chosen based upon a Wilcoxon test
### comparison between objective function values of lambda=0 and
### each lambda>0. The penalized objective function should not
### differ a lot from the unpenalized one. Therefore, the 
### highest lambda that leads to a non-significant difference
### in means is chosen.

  ### Calling SEP.
  funk <- function(a,b){
    x   <- .C("sep",
              as.double(a),
              as.integer(length(a)),
              as.double(b),
              e = integer(length(a)),
              f = double(1),
              PACKAGE="twilight"
              )$f
    return(x)
  }

  ### Try some values.
  lam <- seq(0,0.05,by=0.005)
  n    <- min(1000,length(xin))
  
  boot.p <- matrix(NA,n,50)
  for (i in 1:50){
    boot.p[,i] <- sort(sample(xin,n,replace=TRUE))
  }

  obj0 <- apply(boot.p,2,funk,0)

  if (verbose){
    cat("Wait for some dots. \n")
  }

  res <- NULL
  for (j in 2:length(lam)){
    
    obj.fun <- apply(boot.p,2,funk,lam[j])    
    obj.fun <- obj.fun - obj0
    wilc <- suppressWarnings(wilcox.test(x=obj.fun,mu=0)$p.value)

    ### Look for non-significant differences.
    if (wilc<=0.05){
      res <- lam[j-1]
      break;
    }
        
    if (verbose){cat(".")}
  }
  if (verbose){cat("\n")}

  if (is.null(res)){res <- lam[length(lam)]}
  return(res) 
}
