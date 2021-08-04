plot.twilight <- function(x, which=NULL, grayscale=FALSE, legend=TRUE, ...){
### Plotting function for objects of class "twilight".
### Produces these plots:
###
### "scores": Expected vs. observed test statistics with
###          confidence lines. Points exceeding the confidence lines
###          are highlighted. Reference:
###          Tusher VG, Tibshirani R and Chu G (2001): Significance
###          analysis of mircroarrays applied to the ionizing response,
###          PNAS 98(9), pp. 5116-5121.
###
### "qvalues": q-values vs. number of rejected hypotheses.
###
### "fdr": p-values vs. 1 - local false discovery rate.
###          Bottom ticks are 1% quantiles of p-values.
###          If computed, with bootstrap estimate and bootstrap confidence interval.
###
### "volcano": Volcano plot: Observed scores vs. local false discovery rate.
###          Bottom ticks are 1% quantiles of scores.
###
### "effectsize": Effect size distribution in terms of fold change equivalent scores.
### "table": Tabulate effect size distribution.
###
### Additional input:
### "grayscale": TRUE or FALSE. FALSE produces colored plots.
### "legend":    TRUE or FALSE. Produces legends in "scores", "fdr" and "effectsize".


  funk0 <- function(yin,kol,leg,...){
    if (all(is.nan(yin$result$fdr[1]),na.rm=FALSE)){
      par(mfrow=c(1,2),mar=c(4,4,2,2)+0.1)
      funk1(x,kol,leg,...)
      funk2(x,...)      
    }
    else {
      par(mfrow=c(2,2),mar=c(4,4,2,2)+0.1)
      funk1(x,kol,leg,...)
      funk2(x,...)      
      funk3(x,kol,leg,...)
      funk4(x,...)
    }
  }


  
  funk1 <- function(yin,kol,leg,...){
    if (all(is.nan(yin$result$observed[1]),na.rm=FALSE)){
      stop("The input object must contain observed and expected test scores.\n Choose 'qvalues' or 'fdr' instead.\n")
    }

    maxi <- 2*abs(yin$result$observed[1])
    y <- abs(yin$result$observed[1])
        
    if (kol==TRUE){
      plot(yin$result$expected,yin$result$observed,xlab="Expected score",ylab="Observed score",type="n",xlim=c(-y,y),ylim=c(-y,y),...)
      lines(c(-maxi,maxi),c(-maxi,maxi))
      points(yin$result$expected[as.logical(!yin$result$candidate)],yin$result$observed[as.logical(!yin$result$candidate)])
      points(yin$result$expected[as.logical(yin$result$candidate)],yin$result$observed[as.logical(yin$result$candidate)],col=gray(0.5))
      lines(c(-maxi,maxi),c(-maxi-yin$ci.line,maxi-yin$ci.line),col=gray(0.5))
      lines(c(-maxi,maxi),c(-maxi+yin$ci.line,maxi+yin$ci.line),col=gray(0.5))      
      if (leg==TRUE){
        legend(0,-y*0.8,legend=paste(as.character(100*yin$quant.ci),"% CI",sep=""),lty=1,bty="n",col=gray(0.5))
      }
    }
    
    if (kol==FALSE){
      plot(yin$result$expected,yin$result$observed,xlab="Expected score",ylab="Observed score",type="n",xlim=c(-y,y),ylim=c(-y,y),...)
      lines(c(-maxi,maxi),c(-maxi,maxi),lwd=2)
      points(yin$result$expected[as.logical(!yin$result$candidate)],yin$result$observed[as.logical(!yin$result$candidate)])
      points(yin$result$expected[as.logical(yin$result$candidate)],yin$result$observed[as.logical(yin$result$candidate)],col="red")
      lines(c(-maxi,maxi),c(-maxi-yin$ci.line,maxi-yin$ci.line),col="red")
      lines(c(-maxi,maxi),c(-maxi+yin$ci.line,maxi+yin$ci.line),col="red")      
      if (leg==TRUE){
        legend(0,-y*0.8,legend=paste(as.character(100*yin$quant.ci),"% CI",sep=""),lty=1,bty="n",col="red")
      }
    }
  }





  funk2 <- function(yin,...){
    y        <- unique(yin$result$qvalue)
    hist.num <- hist(yin$result$qvalue,br=c(-1,y),plot=FALSE)$counts

    plot(c(0,y),cumsum(c(0,hist.num)),t="s",xlab="Q-value",ylab="No. of significant outcomes",...)
    lines(c(-10,10),c(0,0),col=gray(0.5))
  }






  funk3 <- function(yin,kol,leg,...){
    if (all(is.nan(yin$result$fdr[1]),na.rm=FALSE)){
      stop("The input object must contain local FDR values.\n Choose 'scores' or 'qvalues' instead or run twilight.\n")
    }

    q.tick <- quantile(yin$result$pvalue,seq(0,1,by=0.01))
    x <- round(seq(1,nrow(yin$result),length=100))
    
    if (kol==TRUE){
      plot(yin$result$pvalue[x],1-yin$result$fdr[x],t="n",xlab="P-value",ylab=expression("1-"~~widehat(fdr)),ylim=c(0,1),...)
      lines(c(-10,10),c(0,0),col=gray(0.5))
      if (is.nan(yin$result$mean.fdr[1])==FALSE){
        lines(yin$result$pvalue[x],1-yin$result$lower.fdr[x],col=gray(0.5),lty=2)
        lines(yin$result$pvalue[x],1-yin$result$upper.fdr[x],col=gray(0.5),lty=2)
        lines(yin$result$pvalue[x],1-yin$result$mean.fdr[x],col=gray(0.5))
        if (leg==TRUE){
          legend(0.6,1,legend=c(expression("1-"~~widehat(fdr)),"Mean",paste(as.character(100*yin$boot.ci),"% CI",sep="")),lty=c(1,1,2),bty="n",col=c("black",gray(0.5),gray(0.5)),lwd=c(2,1,1))
        }
      }
      lines(yin$result$pvalue[x],1-yin$result$fdr[x],lwd=2)
      rug(q.tick)
    }

    if (kol==FALSE){
      plot(yin$result$pvalue[x],1-yin$result$fdr[x],t="n",xlab="P-value",ylab=expression("1-"~~widehat(fdr)),ylim=c(0,1),...)
      lines(c(-10,10),c(0,0),col=gray(0.5))
      if (is.nan(yin$result$mean.fdr[1])==FALSE){
        lines(yin$result$pvalue[x],1-yin$result$lower.fdr[x],col="red",lty=2)
        lines(yin$result$pvalue[x],1-yin$result$upper.fdr[x],col="red",lty=2)
        lines(yin$result$pvalue[x],1-yin$result$mean.fdr[x],col="red")
        if (leg==TRUE){
          legend(0.6,1,legend=c(expression("1-"~~widehat(fdr)),"Mean",paste(as.character(100*yin$boot.ci),"% CI",sep="")),lty=c(1,1,2),bty="n",col=c("black","red","red"),lwd=c(2,1,1))
        }
      }
      lines(yin$result$pvalue[x],1-yin$result$fdr[x],lwd=2)
      rug(q.tick)
    }
  }





  funk4 <- function(yin,...){
    if (all(is.nan(yin$result$fdr[1]),na.rm=FALSE)){
      stop("The input object must contain local FDR values.\n Choose 'scores' or 'qvalues' instead or run twilight.\n")
    }

    q.tick <- quantile(yin$result$observed,seq(0,1,by=0.01))
    maxi   <- 2*max(abs(yin$result$observed))

    plot(yin$result$observed,1-yin$result$fdr,t="n",xlab="Observed test score",ylab=expression("1-"~~widehat(fdr)),ylim=c(0,1),...)
    lines(c(-maxi,maxi),c(0,0),col=gray(0.5))
    points(yin$result$observed,1-yin$result$fdr)
    rug(q.tick)
  }



  funk5 <- function(yin,leg,...){
    if (all(is.nan(yin$effect[[1]]),na.rm=FALSE)){
      stop("The input object must contain effect size frequencies.\n Choose 'scores' or 'qvalues' instead or run twilight.\n")
    }

    check <- unlist(strsplit(yin$call," "))
    if (check[2]!="fc."){
      stop("Effect size distributions can only be estimated with twilight.pval(.,method='fc').\n")
    }

    all <- hist(yin$result$observed,br=yin$effect$breaks,plot=FALSE)

    x.fc   <- yin$result$observed

    mini   <- ceiling(min(x.fc))
    maxi   <- floor(max(x.fc))
    x.tick <- seq(mini,maxi,length=abs(mini)+abs(maxi)+1)

    x.lab  <- (exp(abs(x.tick))-1)*100
    x.lab  <- paste(sign(x.tick)*round(x.lab),"%",sep="")

    plot(all,col=gray(0.7),xaxt="n",main="",xlab="Fold change equivalent score")
    plot(yin$effect,col="black",add=TRUE)

    if (leg==TRUE){
      legend(mean(x.tick[5:6]),max(all$counts),legend=c("Mixture","Alternative"),lty=c(1,1),bty="n",col=c(gray(0.7),"black"),lwd=c(2,2),y.intersp=2)
    }
      
    axis(1,at=x.tick,labels=x.lab)

  }






  funk6 <- function(yin){
    if (all(is.nan(yin$effect[[1]]),na.rm=FALSE)){
      stop("The input object must contain effect size frequencies.\n Choose 'scores' or 'qvalues' instead or run twilight.\n")
    }

    check <- unlist(strsplit(yin$call," "))
    if (check[2]!="fc."){
      stop("Effect size distributions can only be estimated with twilight.pval(.,method='fc').\n")
    }
    
    all <- hist(yin$result$observed,br=yin$effect$breaks,plot=FALSE)
    fce <- round((exp(abs(yin$effect$mids))-1)*100*sign(yin$effect$mids))
    tab <- cbind(yin$effect$mids,all$counts,yin$effect$counts)

    rownames(tab) <- paste(fce,"%",sep="")
    colnames(tab) <- c("LogRatio","Mixture","Alternative")
    
    return(tab)

  }

  


  if (is.null(which)){funk0(x,grayscale,legend,...)}
  else {
    switch(which,
           scores = funk1(x,grayscale,legend,...),
           qvalues = funk2(x,...),
           fdr = funk3(x,grayscale,legend,...),
           volcano = funk4(x,...),
           effectsize = funk5(x,legend,...),
           table = funk6(x)
           )
  }
  
}
