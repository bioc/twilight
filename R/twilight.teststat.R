twilight.teststat <- function(xin,yin,method="fc",paired=FALSE,s0=NULL){
  ### extract data matrix if class(xin) is exprSet
  xin <- twilight.getmatrix(xin)

  ### check dimensions
  if (!is.matrix(xin)){
    stop("First input must be a matrix. \n")
  }

  if (!is.vector(yin)){
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

  if ((method!="pearson")&(method!="spearman")){
  ### translate index vector yin into binary vector with 1 as case and 0 as control samples.
    y <- sort(unique(yin),na.last=TRUE)
    
    if (length(y)!=2){
      stop("Samples must belong to TWO classes. \n")
    }
    
    y1 <- which(yin==y[1])
    y2 <- which(yin==y[2])
    yin[y1] <- 0
    yin[y2] <- 1
    
    if (paired){
      if (sum(yin)!=sum(1-yin)){
        stop("This is a PAIRED twosample test. The sizes of the two classes must be equal. \n")
      }
    }
  }

  ### transform to ranks for Spearman coefficient.
  if (method=="spearman"){
    yin <- rank(yin)
    xin <- t(apply(xin,1,rank))
  }  

  
  ### Z test with s0=0 is a t test.
  if (!is.null(s0)){
    if((s0==0)&(method=="z")){method <- "t"}
  }
  
  if (is.null(s0)){s0 <- 0}


  if ((method!="pearson")&(method!="spearman")){
    funk <- function(a, b, c, d, s) {
      .C(ifelse(paired,"paired","unpaired"), 
         as.integer(a),
         as.integer(sum(a)),
         as.integer(sum(1 - a)),       
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         as.integer(c), 
         as.integer(which(d == 1) - 1),
         as.integer(which(d == 0) - 1),
         as.double(s),
         e = double(nrow(b)), PACKAGE = "twilight")$e
    }
    
    stat.obs <- switch(method,
                       t = funk(yin,xin,1,yin,s0),
                       z = funk(yin,xin,2,yin,s0),
                       fc = funk(yin,xin,3,yin,s0)
                       )
  }
  
  if ((method=="pearson")|(method=="spearman")){
    funk <- function(a,b){
      .C("corsingle",
         as.double(a),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         e=double(nrow(b)),PACKAGE="twilight")$e
    }
    
    stat.obs <- funk(yin,xin)
  }


  return(stat.obs)
}
