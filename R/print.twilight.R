print.twilight <- function(x, ...){

  cat("\n Twilight object with")
  cat("\n    ",nrow(x$result),"transcripts")
  cat("\n     observed and expected test statistics")
  cat("\n     p- and q-values")

  if (!is.nan(x$result$fdr[1])){
    cat("\n     local FDR")    
  }
  if (!is.nan(x$result$mean.fdr[1])){
    cat("\n     bootstrap estimates of local FDR")    
  }

  cat("\n")
  if (is.nan(x$boot.pi0[[1]])){
    cat("\n Estimated percentage of non-induced genes:\n")
    a <- x$pi0
    names(a) <- "pi0"
    print(a)
  }
  if (!is.nan(x$boot.pi0[[1]])){
    cat("\n Bootstrap estimate of percentage of non-induced")
    cat("\n genes with lower and upper ",x$boot.ci*100,"% CI:\n",sep="")
    print(x$boot.pi0)
  }
      
  cat("\n Function call:\n",x$call,"\n")
  if (!is.nan(x$lambda)){
    cat(" Function twilight used lambda =",x$lambda,"\n")
  }
}
