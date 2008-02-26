twilight.getmatrix <- function(xin){
  ## Extracts the data matrix from an expression set.
  if (class(xin)[1]=="ExpressionSet"){
    xout <- exprs(xin)
    rownames(xout) <- featureNames(xin)
    colnames(xout) <- sampleNames(xin)
    result <- xout
  } 
  else {
	result <- xin
  }
  return(result)
}
