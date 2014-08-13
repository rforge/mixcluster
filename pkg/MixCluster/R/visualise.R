compute_w <- function(x, param, k, kind=NULL){
  if (is.null(kind)){
    kind <- rep(1,ncol(x))
    for (j in 1:ncol(x)){
      if (all(x[,j]==ceiling(x[,j]))){
        if (max(x[,j])==1){
          kind[j] <- 3
        }else{
          kind[j] <- 2
        }
      }
    }
  }
  
  y <- x
  if (any(kind==1)){
    for (j in which(kind==1))
      y[,j] <- (x[,j]-param@margins[[k]][j,1])/param@margins[[k]][j,2]
  }
  
  if (any(kind!=1)){
    who <- which(kind!=1)
    for (j in who){
      if (kind[j]==3){
        y[,j] <- genere_y_entier(x[,j],param@margins[[k]][j,1])
      }else if (kind[j]==3){
        y[,j] <- genere_y_binaire(x[,j],param@margins[[k]][j,1]) 
      }
    }
  }
  
  
  sauve <- y
  for (it in 1:500){
    if (any(kind!=1)){
      for (j in who){
        me <-  rowSums(sweep(as.matrix(y[,-j]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
        sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )
        
        
        if (kind[j]==2){
          y[,j] <- genere_y_entier2(param@margins[[k]][j,1], x[,j], me, sig)
        }else if (kind[j]==3){
          y[,j] <- genere_y_binaire2(param@margins[[k]][j,1], x[,j], me, sig) 
        }
      }
    }
    sauve <- sauve + y
  }
  sauve <- as.matrix(sauve /1000)
  return(sauve)
}

compute_cik <- function(obj, class){
  w <- compute_w(obj@data$data, obj@param, class, obj@data$kind)
  res.eigen <- eigen(obj@param@correlations[[class]]) 
  return(w%*%res.eigen$vectors)
}



drawn_cik <- function(obj, class, axe){
  res.eigen <- eigen(obj@param@correlations[[class]]) 
  plot(obj@expect_y[[class]][,axe[1]], obj@expect_y[[class]][,axe[2]],cex.lab=1.5,cex=1,lwd=1.5, col=obj@partition, pch=class+15, 
       xlab=paste("inertia",round(res.eigen$values[axe[1]]/obj@data$ncol*100,1),"%"), ylab=paste("inertia",round(res.eigen$values[axe[2]]/obj@data$ncol*100,1),"%"))
  abline(v=0,lty=2)
  abline(h=0,lty=2)
}

plot_correl <- function(obj, class, axe){
  res.eigen <- eigen(obj@param@correlations[[class]]) 
  plot(NA,xlim=c(-1,1),ylim=c(-1,1),cex.lab=1.5,cex=1,lwd=1.5,
       xlab=paste("inertia",round(res.eigen$values[axe[1]]/obj@data$ncol*100,1),"%"), ylab=paste("inertia",round(res.eigen$values[axe[2]]/obj@data$ncol*100,1),"%"))
  plan <- res.eigen$vectors
  for (j in 1:nrow(plan)){
    arrows(0,0,plan[j,axe[1]],plan[j,axe[2]],lwd=2)
    c1 <- plan[j,axe[1]]-0.05+0.1*(plan[j,axe[1]]>0)
    c2 <- plan[j,axe[2]]-0.05+0.1*(plan[j,axe[2]]>0)
    text(c1,c2,colnames(obj@data$data)[j],cex=1)
  }
  lines(c(0,0),c(-2,2),lty=2)
  lines(c(-2,2),c(0,0),lty=2)
  draw.circle(0, 0, radius = 1,lwd=1) 
}

MixClusvisu <- function(obj, class, axe=c(1,2)){
  op <- par(no.readonly = TRUE)
  
  par(mar=c(4.2,4.2,1,1))
  if (obj@test_expect_y[class]==0){
    cat("  computation of the coordinates\n")
    obj@test_expect_y[class] <- 1
    obj@expect_y[[class]]  <- compute_cik(obj, class)
  }

  cat("  scatterplot\n")
  drawn_cik(obj, class, axe)
  cat("  correlation circle")
  par(mar=c(4.2,3.2,1,1),mfrow=c(1,1),cex.axis=1.5,pty="s")
  plot_correl(obj, class, axe)
  
  
  
  par(op)
  return(obj)
}



plot_continue <- function(object, j){
  main <- colnames(object@data$data)[j]
  hist(object@data$data[,j], freq=FALSE, main=main, cex.main=2.2, plot=TRUE, xlab="", ylab="",
       xlim=c(min(object@data$data[,j]), max(object@data$data[,j])))
  se <- seq(min(object@data$data[,j]),max(object@data$data[,j]),length.out=100)
  tmp <- se*0
  for (k in 1:length(object@param@proportions)){
    tmp <- tmp + object@param@proportions[k]*dnorm(se, object@param@margins[[k]][j,1], object@param@margins[[k]][j,2])
    lines(se, object@param@proportions[k]*dnorm(se, object@param@margins[[k]][j,1], object@param@margins[[k]][j,2]),col=k,lwd=6)
  }
  lines(se,tmp,col="grey",lwd=6,lty=2)
}

plot_entier <- function(object, j){
  main <- colnames(object@data$data)[j]
  se <- min(object@data$data[,j]):max(object@data$data[,j])
  re <- matrix(0,length(object@param@proportions),length(se))
  colnames(re) <- min(se):max(se)
  
  for (k in 1:length(object@param@proportions)){
    re[k,] <- object@param@proportions[k]*dpois(se, object@param@margins[[k]][j,1])
  }
  barplot(re,col=1:length(object@param@proportions),main=main,cex.main=2.2,
          plot=TRUE,xlab="",xlim=c(min(se),max(se)+3))
}


plot_binaire <- function(object, j){
  main <- colnames(object@data$data)[j]
  proba <- matrix(object@param@proportions,length(object@param@proportions),ncol=2)
  for (k in 1:nrow(proba))
    proba[k,] <- proba[k,]*(c(object@param@margins[[k]][j,1],1-object@param@margins[[k]][j,1]))
  colnames(proba) <- c("Prob.0", "Prob.1")
  mp <- barplot(proba, beside = TRUE, axisnames = TRUE, names.arg=names(t),col=1:nrow(proba),
                cex.main=2.5)
  # add unconditional frequencies
  title(main,cex.main=2.2)
  # Add the individual bar labels
  mtext(1, at = mp, text = paste("C",1:nrow(proba)), line = 0, cex = 1)
}

plot_ordinal <- function(object, j){
  cat("TODO")
  plot(NA,xlim=c(0,1),ylim=c(0,1))
}


setMethod(
  f="plot",
  signature=c("MixClusresults"),
  function(x, y, ...){
    if (missing(y)){
      y <- 1: x@data$ncol
    }
    
    if (is.numeric(y)){
      if (max(y)>x@data$ncol)
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
      if (x@data$kind[j]==1){
        plot_continue(x, j)
      }else if (x@data$kind[j]==2){
        plot_entier(x, j)
      }else if (x@data$kind[j]==3){
        plot_binaire(x, j)
      }else if (x@data$kind[j]==4){
        plot_ordinal(x, j)
      }
    }
    close.screen(all.screens = TRUE)
    # restore plotting parameters
    par(op)
    
  }
  
  
)
