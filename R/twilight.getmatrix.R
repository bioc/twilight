twilight.getmatrix <- function(xin){
### Extracts the data matrix from an expression set.

  if (class(xin)[1]=="exprSet"){
    xout <- exprs(xin)
    rownames(xout) <- geneNames(xin)
    colnames(xout) <- sampleNames(xin)
    return(xout)
  } 
  if (class(xin)[1]=="ExpressionSet"){
    xout <- exprs(xin)
    rownames(xout) <- featureNames(xin)
    colnames(xout) <- sampleNames(xin)
    return(xout)
  } 
  if ((class(xin)[1]!="ExpressionSet")&(class(xin)[1]!="exprSet")){return(xin)}
}
