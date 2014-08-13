GCMMpriors <- function(x,kind){
  priors <- list();
  if (any(kind==1)){
    for (j in which(kind==1)){
      priors[[j]] <- list()
      priors[[j]]$bo<- mean(x[,j])
      priors[[j]]$No <- 2.6/(max(x[,j])-min(x[,j]))
      priors[[j]]$co <- 1.28
      priors[[j]]$Co <- 0.36*var(x[,j])
    }
  }
  return(priors)  
}

compute_proba <- function(x,param,kind){
  proba <- matrix(param@proportions,nrow(x),length(param@proportions),byrow=TRUE)
  for (k in 1:ncol(proba)){
    if (any(kind==1)){
      y <- matrix(0,nrow(x),sum(kind==1))
      cp <- 1
      for (j in which(kind==1)){
        y[,cp] <- (x[,j]-param@margins[[k]][j,1])/param@margins[[k]][j,2]
        cp <- cp + 1
      }
      if (sum(kind==1)>1){
        proba[,k] <- proba[,k] * dmvnorm(y,rep(0,sum(kind==1)),sigma=param@correlations[[k]][which(kind==1),which(kind==1)])/prod(param@margins[[k]][which(kind==1),2])
      }else{
        proba[,k] <- proba[,k] * dnorm(y,0,1)/param@margins[[k]][which(kind==1),2]  
      }
    }
    if (any(kind>1)){
      for (i in 1:nrow(x)){
        binf <- rep(0,sum(kind!=1))
        bsup <- binf
        cp <- 1
        for (j in which(kind!=1)){
          if (kind[j]==2){
            binf[cp] <- qnorm(ppois(x[i,j]-1,param@margins[[k]][j,1]))
            bsup[cp] <- qnorm(ppois(x[i,j],param@margins[[k]][j,1]))
          }else if (kind[j]==3){
            if (x[i,j]==0){
              binf[cp] <- -Inf
              bsup[cp] <- qnorm(param@margins[[k]][j,1])
            }else{
              binf[cp] <- qnorm(param@margins[[k]][j,1])
              bsup[cp] <- Inf
            }
          }
          cp <- cp + 1
        }
        if (all(binf<bsup)){
          if (any(kind==1)){
            me <- (param@correlations[[k]][which(kind!=1),which(kind==1)] %*% solve(param@correlations[[k]][which(kind==1),which(kind==1)]) ) %*% y[i,]
            sig <- param@correlations[[k]][which(kind!=1),which(kind!=1)] - param@correlations[[k]][which(kind!=1),which(kind==1)] %*% solve(param@correlations[[k]][which(kind==1),which(kind==1)])  %*%param@correlations[[k]][which(kind==1),which(kind!=1)]
          }else{
            me <- rep(0,sum(kind>1))
            sig <- param@correlations[[k]]
          }
          proba[i,k] <- proba[i,k] * pmvnorm(binf,bsup,mean=as.numeric(me),sigma=sig)
        }else{
          proba[i,k] <- 0
        }
        
      }
    }    
  }
  return(proba)
}

compute_criteria <- function(x,param,kind,hetero=TRUE, cim=FALSE, uncoeff=FALSE){
  proba <- as.matrix(compute_proba(x,param,kind))
  tik <- proba/rowSums(proba)
  z <- rep(0,nrow(tik))
  icl <- 0
  for (k in 1:ncol(tik)){
    z[which(rowSums(sweep(tik,1,tik[,k],"-")<=0)==ncol(tik))] <- k
  }
  for (k in 1:ncol(tik)){
    icl <- icl + sum(log(proba[which(z==k),k]))
  }
  
  loglike <- sum(log(rowSums(proba)))
  g <- ncol(proba)
  d <- ncol(x)
  if (hetero==TRUE){
    if (uncoeff==FALSE){
      nbparam <- (g-1)+g*(d*(d-1)/2 + d)+g*sum(kind==1)
    }else{
      nbparam <- (g-1)+g*(1 + d)+g*sum(kind==1)
    }
    
    
  }else{
    if (cim==FALSE){
      if (uncoeff==FALSE){
        nbparam <- (g-1)+ d*(d-1)/2 + d*g +g*sum(kind==1)
      }else{
        nbparam <- (g-1)+ 1 + d*g +g*sum(kind==1)
      }
    }else{
      nbparam <- (g-1) + d*g +g*sum(kind==1)
      
    }
  }
  
  bic <- loglike - nbparam*0.5*log(nrow(x))
  icl <- icl - nbparam*0.5*log(nrow(x))
  return(list(tik=tik, loglike=loglike, bic=bic, icl=icl, nbparam=nbparam, z=z))
}


addparam <- function(save_param,param){
  save_param@proportions <- save_param@proportions + param@proportions
  for (k in 1:length(save_param@proportions)){
    save_param@margins[[k]] <- save_param@margins[[k]] + param@margins[[k]]
    save_param@correlations[[k]] <- save_param@correlations[[k]] + param@correlations[[k]]
  }
  return(save_param)
}

normsave <- function(param,itermax){
  param@proportions <- param@proportions/(itermax+1)
  for (k in 1:length(param@proportions)){
    param@margins[[k]] <- param@margins[[k]]/(itermax+1)
    param@correlations[[k]] <- param@correlations[[k]]/(itermax+1)
  }
  return(param)
}