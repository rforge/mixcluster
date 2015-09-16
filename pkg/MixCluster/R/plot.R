
plot_continue <- function(object, j){
  main <- colnames(object@data@x)[j]
  e <- range(object@data@x[object@data@o[[j]], j])
  hist(object@data@x[,j], freq=FALSE, main=main, cex.main=2.2, plot=TRUE, xlab="", ylab="",
       xlim=e)
  se <- seq(e[1], e[2],length.out=100)
  tmp <- se*0
  for (k in 1:length(object@param@pi)){
    tmp <- tmp + object@param@pi[k]*dnorm(se, object@param@beta[[k]][[j]]@mu, object@param@beta[[k]][[j]]@sd)
    lines(se, object@param@pi[k]*dnorm(se, object@param@beta[[k]][[j]]@mu, object@param@beta[[k]][[j]]@sd), col=k, lwd=6)
  }
  lines(se,tmp,col="grey",lwd=6,lty=2)
}

plot_entier <- function(object, j){
  main <- colnames(object@data@x)[j]
  e <- range(object@data@x[object@data@o[[j]], j])+c(-10,10)
  se <- e[1]: e[2]
  re <- matrix(0,length(object@param@pi),length(se))
  colnames(re) <- min(se):max(se)
  
  for (k in 1:length(object@param@pi))
    re[k,] <- object@param@pi[k]*dpois(se, object@param@beta[[k]][[j]]@beta)
  
  barplot(re,col=1:length(object@param@pi),main=main,cex.main=2.2,
          plot=TRUE,xlab="",xlim=c(min(se),max(se)+3))
}

plot_ordinal <- function(object, j){
  main <- colnames(object@data@x)[j]
  proba <- matrix(object@param@pi,length(object@param@pi),ncol=length(unique(object@data@x[object@data@o[[j]],j])))
  for (k in 1:nrow(proba))
    proba[k,] <- proba[k,]*object@param@beta[[k]][[j]]@beta
  colnames(proba) <- paste("Prob.", 0:(ncol(proba)-1),sep="")
  mp <- barplot(proba, beside = TRUE, axisnames = TRUE, names.arg=names(t),col=1:nrow(proba),
                cex.main=2.5)
  # add unconditional frequencies
  title(main,cex.main=2.2)
  # Add the individual bar labels
  mtext(1, at = mp, text = paste("C",1:nrow(proba)), line = 0, cex = 0.6)
}


setMethod(
  f="plot",
  signature=c("MixClusResults"),
  function(x, y){
    if (missing(y)){
      y <- 1: x@data@e
    }
    
    if (is.numeric(y)){
      if (max(y)>x@data@e)
        stop("y mismatchs the data frame dimension")
    }else{
      stop("y has to be a numeric vector indicating the variables to be plotted")
    }
    # get old par 
    op <- par(no.readonly = TRUE) 
    
    par(mar = rep(2.5,4),cex = .75)
    if (length(y)==1){
      split.screen(c(1,1))
    }else if (length(y)<=3){
      split.screen(c(1,length(y)))
    }else if (length(y)==4){
      split.screen(c(2,2))
    }else if (length(y)<13){
      split.screen(c(length(y)%/%3 + 1*(length(y)!= (3*(length(y)%/%3))),3))
    }else{
      split.screen(c(length(y)%/%4 + 1*(length(y)!= (4*(length(y)%/%4))),4))
    }
    cp <- 0
    for ( j in y){
      cp <- cp+1
      screen(cp)
      if (x@data@kind[j]==1){
        plot_continue(x, j)
      }else if (x@data@kind[j]==2){
        plot_entier(x, j)
      }else if (x@data@kind[j]==3){
        plot_ordinal(x, j)
      }
    }
    close.screen(all.screens = TRUE)
    # restore plotting parameters
    par(op)
    
  }
  
  
)
