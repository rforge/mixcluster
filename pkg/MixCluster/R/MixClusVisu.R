MixClusUpdateForVisu <- function(output){
  for (k in 1:output@model$model@g){
    output@data@condexpec[[k]] <- matrix(0, output@data@n, output@data@e)
    sup <- matrix(0,output@data@n,output@data@e)
    inf <- matrix(0,output@data@n,output@data@e)
    for (j in 1:output@data@e){
      bound <- cbind(rep(-Inf, output@data@n), rep(Inf, output@data@n))
      bound[output@data@o[[j]], ] <- findbounds(output@data@x[output@data@o[[j]], j], output@param@beta[[k]][[j]])
      if (output@data@kind[j]==1){
        bound[,1] <- bound[,1]-10**(-6) 
      }
      sup[,j] <- bound[,2]
      if (any(sup[,j]==-Inf)) sup[which(sup[,j]==-Inf),j] <- -10**(6)
      inf[,j] <- bound[,1]
      if (any(inf[,j]==Inf)) inf[which(inf[,j]==Inf),j] <- 10**(6)
    }
    for (i in 1:output@data@n){
      print(i)
      output@data@condexpec[[k]][i,] <- mtmvnorm(mean=rep(0,output@data@e), sigma=output@param@correl[[k]], lower=inf[i,], upper=sup[i,], doComputeVariance = FALSE)$tmean
    }
  
    
    
  }
       return(output)
}

MixClusVisu <- function(obj, class, axe=c(1,2), figure=c("scatter","circle"),...){
  op <- par(no.readonly = TRUE)
  res.eigen <- eigen(obj@param@correl[[class]]) 
  w <- obj@data@condexpec[[class]]%*%res.eigen$vectors
  
  
  if( any(figure=="scatter")){
    
    par(mar=c(4.2,4.2,1,1))
       cat("  scatterplot\n")
    plot(w[,axe[1]], w[,axe[2]],cex.lab=1.5,cex=1,lwd=1.5, col=obj@data@partition, pch=obj@data@partition+15, 
         xlab=paste("inertia",round(res.eigen$values[axe[1]]/obj@data@e*100,1),"%"),
         ylab=paste("inertia",round(res.eigen$values[axe[2]]/obj@data@e*100,1),"%"),...)
    abline(v=0,lty=2)
    abline(h=0,lty=2)
  }


  if (any(figure=="circle")){
    cat("  correlation circle")
    par(mar=c(4.2,3.2,1,1),mfrow=c(1,1),cex.axis=1.5,pty="s")
    
    plot(NA,xlim=c(-1,1),ylim=c(-1,1),cex.lab=1.5,cex=1,lwd=1.5,
         xlab=paste("inertia",round(res.eigen$values[axe[1]]/obj@data@e*100,1),"%"), 
         ylab=paste("inertia",round(res.eigen$values[axe[2]]/obj@data@e*100,1),"%"))
    plan <- res.eigen$vectors
    for (j in 1:nrow(plan)){
      arrows(0,0,plan[j,axe[1]],plan[j,axe[2]],lwd=2)
      c1 <- plan[j,axe[1]]-0.05+0.1*(plan[j,axe[1]]>0)
      c2 <- plan[j,axe[2]]-0.05+0.1*(plan[j,axe[2]]>0)
      text(c1,c2,colnames(obj@data@x)[j],cex=1)
    }
    lines(c(0,0),c(-2,2),lty=2)
    lines(c(-2,2),c(0,0),lty=2)
    draw.circle(0, 0, radius = 1,lwd=1) 
  }


  
  par(op)
  return(invisible())
}


