twilight.permute.pair <- function(v,m,bal=TRUE){
### do "m" permutations of the vector "v" containing paired group labels.

  x <- sort(unique(v))
  
  if (bal==FALSE){
    M <- matrix(NA,m,length(v))
    M[1,] <- v
    for (i in 2:m){
      r <- runif(length(v)/2)
      s <- numeric(length(v))
      r <- round(r)
      s[v==x[1]] <- r
      s[v==x[2]] <- 1-r
      M[i,] <- s
    }
  }

  if (bal==TRUE){
    M <- matrix(NA,m,length(v))
    for (i in 1:m){
      n <- ceiling(length(v)/4)
      r <- runif(n)
      r <- round(r)
      r <- c(r,(1-r)[1:(length(v)/2-n)])
      s <- numeric(length(v))
      s[v==x[1]] <- r
      s[v==x[2]] <- 1-r      
      M[i,] <- s
    }
  }

  return(M)
}
