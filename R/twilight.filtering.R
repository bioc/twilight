twilight.filtering <- function(xin,yin,method="fc",paired=FALSE,s0=0,verbose=TRUE,num.perm=1000,num.take=50){


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
    
    if (paired==TRUE){
      if (sum(yin)!=sum(1-yin)){
        stop("This is a PAIRED twosample test. \n")
      }
    }
  }

  ### transform to ranks for Spearman coefficient.
  if (method=="spearman"){
    yin <- rank(yin)
    xin <- t(apply(xin,1,rank))
  }  


  ### run length check for iteration
  if ((num.perm/num.take > 50)&(verbose)){
    cat("The filtering might take some time. Consider to increase the stepsize 'num.take'. \n")
  }


  ### translate methods into numerical values
  if (method=="t"){meth=1}
  if (method=="z"){meth=2}
  if (method=="fc"){meth=3}
  if (method=="pearson"){meth=4}
  if (method=="spearman"){meth=4}

  
  ### remove non-unique permutations
  simplify <- function(a){    
    comp <- function(x,y,Z){
      funk <- function(i,a,b,C){(identical(C[a[i],],C[b[i],]))|(identical(C[a[i],],1-C[b[i],]))}
      sapply(seq(along=x),funk,x,y,Z)
    }    
    x <- outer(1:nrow(a),1:nrow(a),comp,a)
    diag(x) <- FALSE
    i <- 1
    while (sum(x[upper.tri(x)])>0){
      a <- a[!x[i,],]
      x <- x[!x[i,],!x[i,]]
      i <- i+1
    }    
    return(a)
  }

  
  ### compute test statistics, pooled p-values and KS statistics at once
  if (paired){
    perm <- function(a,x,y){
      b <- twilight.permute.pair(a,max(y,100),bal=FALSE)
      b <- rbind(b,x)
      b <- simplify(b)
      b <- b[-1,]
      return(b)
    }
    
    funk <- function(a,b,c,orig,s,verbose){
      if (verbose){cat(".")}
      .C("pairedKSTEST",
         as.integer(t(a)),
         as.integer(nrow(a)),
         as.integer(sum(a[1,])),
         as.integer(sum(1-a[1,])),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         as.integer(c),
         as.integer(which(orig==1)-1),
         as.integer(which(orig==0)-1),
         as.double(s),
         f=double(nrow(a)),PACKAGE="twilight"
         )$f
    }
  }
  if (!paired){
    perm <- function(a,x,y){
      b <- twilight.permute.unpair(a,max(y,500),bal=FALSE)
      b <- rbind(b,x)
      b <- simplify(b)
      b <- b[-1,]
      return(b)
    }

    funk <- function(a,b,c,orig,s,verbose){
      if (verbose){cat(".")}
      .C("unpairedKSTEST",
         as.integer(t(a)),
         as.integer(nrow(a)),
         as.integer(sum(a[1,])),
         as.integer(sum(1-a[1,])),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         as.integer(c),
         as.double(s),
         f=double(nrow(a)),PACKAGE="twilight"
         )$f
    }
  }
  if (meth==4){
    perm <- function(a,x,y){
      b <- twilight.permute.unpair(a,max(y,500),bal=FALSE)
      b <- rbind(b,x)
      b <- simplify(b)
      b <- b[-1,]
      return(b)
    }

    funk <- function(a,b,c,orig,s,verbose){
      if (verbose){cat(".")}
      .C("correlationKSTEST",
         as.double(t(a)),
         as.integer(nrow(a)),
         as.double(t(b)),
         as.integer(nrow(b)),
         as.integer(ncol(b)),
         f=double(nrow(a)),PACKAGE="twilight"
         )$f
    }
  }


  ### Expected frequencies of Hamming distances
  hamming.exp <- function(x,paired=TRUE){
    if (!paired){
      ham <- seq(0,min(sum(x),sum(1-x)))
      y <- choose(sum(x),ham)*choose(sum(1-x),ham)
    }
    if (paired){
      ham <- seq(0,floor(sum(x)/2))
      y <- choose(sum(x),ham)
      if (sum(x)%%2==0){
        y[length(y)] <- y[length(y)]/2
      }
    }
    y <- y/sum(y)
    names(y) <- ham*2
    return(y)  
  }

  ### Observed table of Hamming distances
  hamming.obs <- function(a,x,paired=TRUE){
    hamming.dist <- function(w,v){ # contributed by F. Markowetz
      u <- table(w,v)
      u <- u[1,2]+u[2,1]
      return(u)
    }

    if (!paired){
      ham <- seq(0,min(sum(x),sum(1-x)))*2
    }
    if (paired){
      ham <- seq(0,floor(sum(x)/2))*2
    }
    d <- apply(a,1,hamming.dist,v=x)
    d <- table(d)
    y <- numeric(length(ham))
    names(y) <- ham
    y[names(d)] <- d
    return(y)
  }
  

  

  
  ### start the filtering

  if (verbose){cat(paste("Filtering: Wait for",ceiling(num.perm/num.take),"to",round(num.perm/num.take)+10,"dots "))}
  id.perm <- perm(yin,NULL,num.take)
  ks <- funk(id.perm,xin,meth,id,s0,verbose)
  a <- sort(ks,index.return=TRUE)  
  if (num.take==1){id.opt <- id.perm[a$ix[1:2],]}
  if (num.take>1){id.opt <- id.perm[a$ix[1:min(num.take,length(ks))],]}

  j <- 2
  
  repeat {
    id.perm <- perm(yin,id.opt,num.take)
    ks <- funk(id.perm,xin,meth,id,s0,verbose)
    a <- sort(ks,index.return=TRUE)  
    b <- a$ix[1:min(j*num.take,length(ks))]
    id.opt <- id.perm[b,]
    ks <- ks[b]
    j <- j+1
    
    if (nrow(id.opt)>=num.perm){break}
    if (j>ceiling(num.perm/num.take)+10){
      warning("Filtering stopped. Number of filtered permutations will be lower than intended.")
      break
    }
  }
  if (verbose){cat("done\n")}

  ### first row contains original labeling but give back exactly B permutations.
  yperm <- rbind(yin,id.opt[-1,])

  ### do a goodness-of-fit chi^2-test to see if at least the distribution of Hamming distances of the filtered permutations to the original labeling is random.
  #a <- hamming.exp(yin,paired=paired)
  #b <- hamming.obs(id.opt,yin,paired=paired)

  ### concatenate small entries at both ends of the distributions.
  #a <- round(a*sum(b))

  #for (i in 1:2){
  #  x <- which(a>=10)[1]-1
  #  a[x] <- sum(a[1:x])
  #  a <- a[-(1:(x-1))]
  #  a <- rev(a)
  #  b[x] <- sum(b[1:x])
  #  b <- b[-(1:(x-1))]
  #  b <- rev(b)
  #}
                   
  #test <- suppressWarnings(chisq.test(x=b,p=a,rescale.p=TRUE)$p.value)
  
  return(yperm)
}


