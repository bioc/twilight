twilight.permute.unpair <- function(v,m,bal=TRUE){
### do "m" permutations of the vector "v" containing unpaired group labels.

  x <- sort(unique(v))
  
  if (bal==FALSE){
    M <- matrix(NA,m,length(v))
    M[1,] <- v
    for (i in 2:m){
      M[i,] <- sample(v)
    }
  }

  if (bal==TRUE){
    M  <- matrix(NA,m,length(v))

    for (i in 1:(m%/%2)){
      v1 <- sample(which(v==x[1]))
      v2 <- sample(which(v==x[2]))

      n1 <- length(v1)/2
      y1 <- c( rep(x[1],floor(n1))  , rep(x[2],ceiling(n1)) )
      y2 <- c( rep(x[1],ceiling(n1)), rep(x[2],floor(length(v)-n1-ceiling(n1))) )
      y  <- c(y1,y2)
      z  <- c(v1,v2)
      z  <- sort(z,index.return=TRUE)

      M[i,] <- y[z$ix]
    }

    for (i in (m%/%2+1):m){
      v1 <- sample(which(v==x[1]))
      v2 <- sample(which(v==x[2]))

      n1 <- length(v2)/2
      y1 <- c( rep(x[2],floor(n1))  , rep(x[1],ceiling(n1)) )
      y2 <- c( rep(x[2],ceiling(n1)), rep(x[1],floor(length(v)-n1-ceiling(n1))) )
      y  <- c(y1,y2)
      z  <- c(v1,v2)
      z  <- sort(z,index.return=TRUE)

      M[i,] <- y[z$ix]
    }
    
    M[1,] <- v
    rownames(M) <- NULL
  }

  return(M)
}




