twilight.combi <- function(xin,pin,bin){
### Function gives all possible permutations of binary vector
### "xin".
###
### INPUT
###
### "xin": Binary vector of class labels.
### "pin": TRUE or FALSE. Are the samples paired or not?
### "bin": TRUE or FALSE. Should the permutations be balanced or 
###        unbalanced?
###
### OUTPUT
###
### "xout": Matrix of permuted class label vector. Each row contains
###         one permutation. If the sample size exceeds certain values,
###         "NULL" is returned.

  vec.runs <- function(x){
    y <- 1
    for (i in 2:length(x)){
      if (x[i]!=x[i-1]){
        y <- c(y,i)
      }
    }
    return(y)
  }
  
  mat.runs <- function(x){
    y <- 1
    for (i in 2:dim(x)[1]){
      if (sum(x[i,]!=x[i-1,])>0){
        y <- c(y,i)
      }
    }
    return(y)
  }
  
  pos <- function(x){
    if (x<0){x <- 0}
    return(x)
  }
  
  unpair.unbal <- function(xin){
    
    n1 <- sum(xin)
    n0 <- sum(1-xin)
    n  <- n0+n1
    m  <- choose(n,n0)
    
    A     <- matrix(NA,m,n)
    b     <- c(m*n0/n,m*n1/n)
    A[,1] <- rep(c(0,1),b)
    
    for (i in 2:n){
      r <- vec.runs(A[,1])
      if (i>2){r <- mat.runs(A[,1:(i-1)])}
      s1 <- rowSums(A[r,],na.rm=TRUE)			
      s0 <- rowSums((1-A[r,]),na.rm=TRUE)			
      
      b.new <- numeric()
      for (j in 1:length(r)){
        b.new <- c(b.new,b[j]*pos(n0-s0[j])/(n-i+1),b[j]*pos(n1-s1[j])/(n-i+1))
      }
      
      A[,i] <- rep(rep(c(0,1),length(b.new)/2),b.new)
      b <- b.new[b.new>0]
    }
    
    
    return(A)
  }
  
  unpair.bal <- function(xin){
    
    n1 <- sum(xin)
    n0 <- sum(1-xin)
    n  <- n0+n1
    
    if (n0<=n1){
      m.ceil  <- choose(n0,ceiling(n0/2))*choose(n1,floor(n0/2))      
      m.floor <- choose(n0,floor(n0/2))*choose(n1,ceiling(n0/2))
      
      if (n0%%2==0){A <- matrix(NA,m.ceil,n)}
      if (n0%%2!=0){A <- matrix(NA,m.ceil+m.floor,n)}
      
      A1 <- unpair.unbal(rep(c(0,1),c(ceiling(n0/2),floor(n0/2))))
      A2 <- unpair.unbal(rep(c(0,1),c(floor(n0/2),n1-floor(n0/2))))
      
      b1 <- rep(1:dim(A1)[1],rep(dim(A2)[1],m.ceil/dim(A2)[1]))
      b2 <- rep(1:dim(A2)[1],dim(A1)[1])
      
      for (i in 1:m.ceil){
        A[i,] <- c(A1[b1[i],],A2[b2[i],])
      }
      
      if (n0%%2!=0){
        A1 <- unpair.unbal(rep(c(0,1),c(floor(n0/2),ceiling(n0/2))))
        A2 <- unpair.unbal(rep(c(0,1),c(ceiling(n0/2),n1-ceiling(n0/2))))
        
        b1 <- rep(1:dim(A1)[1],rep(dim(A2)[1],m.floor/dim(A2)[1]))
        b2 <- rep(1:dim(A2)[1],dim(A1)[1])
        
        for (i in 1:m.floor){
          A[m.ceil+i,] <- c(A1[b1[i],],A2[b2[i],])
        }
      }
    }
        
    if (n1<n0){
      m.ceil  <- choose(n1,ceiling(n1/2))*choose(n0,floor(n1/2))
      m.floor <- choose(n1,floor(n1/2))*choose(n0,ceiling(n1/2))
      
      if (n1%%2==0){A <- matrix(NA,m.ceil,n)}
      if (n1%%2!=0){A <- matrix(NA,m.ceil+m.floor,n)}
      
      A1 <- unpair.unbal(rep(c(0,1),c(n0-floor(n1/2),floor(n1/2))))
      A2 <- unpair.unbal(rep(c(0,1),c(floor(n1/2),ceiling(n1/2))))
      
      b1 <- rep(1:dim(A1)[1],rep(dim(A2)[1],m.ceil/dim(A2)[1]))
      b2 <- rep(1:dim(A2)[1],dim(A1)[1])
      
      for (i in 1:m.ceil){
        A[i,] <- c(A1[b1[i],],A2[b2[i],])
      }
      
      if (n1%%2!=0){
        A1 <- unpair.unbal(rep(c(0,1),c(n0-ceiling(n1/2),ceiling(n1/2))))
        A2 <- unpair.unbal(rep(c(0,1),c(ceiling(n1/2),floor(n1/2))))
        
        b1 <- rep(1:dim(A1)[1],rep(dim(A2)[1],m.floor/dim(A2)[1]))
        b2 <- rep(1:dim(A2)[1],dim(A1)[1])
        
        for (i in 1:m.floor){
          A[m.ceil+i,] <- c(A1[b1[i],],A2[b2[i],])
        }
      }
    }
    
    return(A)
  }
  
  pair.unbal <- function(xin){
    
    n <- sum(xin)
    m <- 2^(n-1)
    A <- matrix(NA,m,2*n)
    
    b2 <- numeric()
    for (i in 1:floor(n/2)){
      b1 <- rep(c(0,1),c(n-i,i))
      b2 <- rbind(b2,unpair.unbal(b1))    
    }
    
    A[1,] <- sort(xin)
    for (i in 2:m){
      A[i,] <- c(b2[i-1,],1-b2[i-1,])
    }
    
    return(A)
  }
  
  pair.bal <- function(xin){
    
    n <- sum(xin)
    
    m.floor <- choose(n,floor(n/2))/2
    m.ceil  <- choose(n,ceiling(n/2))/2
    
    if (n%%2==0){m <- m.floor}
    if (n%%2!=0){m <- m.floor+m.ceil}
    
    A <- matrix(NA,m,2*n)
    
    b1 <- rep(c(0,1),c(n-floor(n/2),floor(n/2)))
    b2 <- unpair.unbal(b1)
    b2 <- b2[1:(dim(b2)[1]/2),]
    
    if (n%%2!=0){
      b1 <- rep(c(0,1),c(n-ceiling(n/2),ceiling(n/2)))
      b2 <- rbind(b2,unpair.unbal(b1))
    }
    
    for (i in 1:m){
      A[i,] <- c(b2[i,],1-b2[i,])
    }
    
    return(A)
  }



  n1 <- sum(xin)
  n0 <- sum(1-xin)
  n  <- n1+n0

  if (n>170){
    xout <- NULL
  }

  if (n<=170){
    if (pin==FALSE){
      if (bin==FALSE){
        xout <- NULL
        m <- choose(n,n0)
        if (m<=10000){
          xout <- unpair.unbal(xin)
          
          if (n0==n1){
            xout <- xout[-dim(xout)[1],]
          }
        }
      }
      if (bin==TRUE){
        xout <- NULL
        
        n.min <- min(n0,n1)
        n.max <- max(n0,n1)
        m <- choose(n.min,floor(n.min/2))*choose(n.max,floor(n.min/2))
        if (n.min%%2!=0){m <- m + choose(n.min,ceiling(n.min/2))*choose(n.max,ceiling(n.min/2))}
        
        if (m<=10000){
          xout <- unpair.bal(xin)
        }
      }
    }
    
    if (pin==TRUE){

      if (n0!=n1){
        stop("The input vector is not paired.")
      }
      
      n <- n/2
      
      if (bin==FALSE){
        xout <- NULL
        m <- 2^(n-1)
        if (m<=10000){
          xout <- pair.unbal(xin)
        }
      }
      if (bin==TRUE){
        xout <- NULL
        m.floor <- choose(n,floor(n/2))/2
        m.ceil  <- choose(n,ceiling(n/2))/2
        
        m <- m.floor
        if (n%%2!=0){m <- m+m.ceil}
        if (m<=10000){
          xout <- pair.bal(xin)
        }
      }
    }
  }
  
  ix <- rank(xin,ties.method="first")
  xout <- xout[,ix]

  if (bin==TRUE){
    xout <- rbind(xin,xout)
    rownames(xout) <- NULL
  }
  
  return(xout)
}
