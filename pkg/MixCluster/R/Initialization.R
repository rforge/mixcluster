GCMMparam_init_EM_loc_indpt <- function(g,x,kind){
  pi <- as.numeric(rdirichlet(1,rep(2,g)))
  margins <- list()
  gamma <- list()
  for (k in 1:g){
    margins[[k]] <- matrix(0,ncol(x),2)
    for (j in 1:ncol(x)){
      if (kind[j]==1){
        margins[[k]][j,] <- c(x[sample(1:nrow(x),1),j],sd(x[,j]))
      }else if (kind[j]==2){
        margins[[k]][j,] <- c(x[sample(1:nrow(x),1),j]+0.2,0)
      }else if (kind[j]==3){
        margins[[k]][j,] <- c(rbeta(1,1/2,1/2),0)
      }
    }
    gamma[[k]] <- diag(1,ncol(x),ncol(x))
  }
  return(GCMMparam(pi,margins,gamma))
}

GCMMparam_init <- function(g,x,kind){
  base <- XEM_loc_indpt(x,g,kind)
  for (u in 1:1){
    cand <- XEM_loc_indpt(x,g,kind)
    if (cand$loglike>base$loglike)
      base <- cand
  }
  return(base$param)
}

GCMMlatent_init <- function(x,g,param,kind){
  latent <- matrix(0,nrow(x),ncol(x)+1)
  proba <- compute_proba_cond(x, g, param, kind)
  proba <- as.matrix(proba/rowSums(proba))
  proba <- t(apply(proba,1,cumsum))
  tmp <- runif(nrow(x))
  if (g>1){
    latent[,1] <- rowSums(sweep(proba,1,tmp,"<="))+1
  }else{
    latent[,1] <- rep(1,nrow(x))
  }

  for (k in unique(latent[,1])){
    for (j in 1:ncol(x)){
      latent[which(latent[,1]==k),j+1] <- sample_yj_indpt(x[which(latent[,1]==k),j],kind[j],param@margins[[k]][j,])
    }
  }
  return(latent)
}


compute_proba_cond <- function(x,g,param,kind){
  proba <- matrix(param@proportions,nrow(x),g,byrow=TRUE)
  for (j in 1:ncol(x)){
    if (kind[j]==1){
      for (k in 1:g){
        proba[,k] <- proba[,k] * dnorm(x[,j],param@margins[[k]][j,1],param@margins[[k]][j,2])
      }
    }else if (kind[j]==2){
      for (k in 1:g){
        proba[,k] <- proba[,k] * dpois(x[,j],param@margins[[k]][j,1])
      }
    }else if (kind[j]==3){
      for (k in 1:g){
        proba[,k] <- proba[,k] * ((x[,j]==0)*param@margins[[k]][j,1]+(x[,j]==1)*(1-param@margins[[k]][j,1]))
      }
    }
  }
  return(proba)
}

XEM_loc_indpt <- function(x,g,kind){
  param <- GCMMparam_init_EM_loc_indpt(g,x,kind)
  proba <- compute_proba_cond(x, g, param, kind)
  loglike <- sum(log(rowSums(proba)))
  prec <- -Inf
  loglike
  while ((loglike-prec)>10^(-2)){
    proba <- proba/rowSums(proba)
    
    param@proportions <- colSums(proba)/nrow(x)
    new_proba <- matrix(param@proportions,nrow(x),g,byrow=TRUE)
    for (k in 1:g){
      for (j in 1:ncol(x)){        
        if (kind[j]==1){
          param@margins[[k]][j,1] <- sum(proba[,k]*x[,j])/sum(proba[,k])
          param@margins[[k]][j,2] <- sqrt(sum(proba[,k]*((x[,j]-param@margins[[k]][j,1])**2))/sum(proba[,k]))
          new_proba[,k] <- new_proba[,k]  * dnorm(x[,j],param@margins[[k]][j,1],param@margins[[k]][j,2])
        }else if (kind[j]==2){
          param@margins[[k]][j,1] <- sum(proba[,k]*x[,j])/sum(proba[,k])
          new_proba[,k] <- new_proba[,k]  * dpois(x[,j],param@margins[[k]][j,1])
        }else if (kind[j]==3){
          param@margins[[k]][j,1] <- 1-sum(proba[,k]*x[,j])/sum(proba[,k])
          new_proba[,k] <- new_proba[,k]  * ((x[,j]==0)*param@margins[[k]][j,1]+(x[,j]==1)*(1-param@margins[[k]][j,1]))
        }
      }
    }
    proba <- new_proba
    #proba <- compute_proba_cond(x, g, param, kind)
    prec <- loglike
    loglike <- sum(log(rowSums(proba)))
  
    if(is.na(loglike)){
      
      param <- GCMMparam_init_EM_loc_indpt(g,x,kind)
      proba <- compute_proba_cond(x, g, param, kind)
      loglike <- sum(log(rowSums(proba)))
    }
    if (prec>loglike){print(c(prec,loglike))}
  }
  
  for (j in 1:ncol(x)){
    for (k in 1:g){
      if (kind[j]==3){
        if (param@margins[[k]][j,1]>0.95)
          param@margins[[k]][j,1] <- 0.95
        if (param@margins[[k]][j,1]<0.05)
          param@margins[[k]][j,1] <- 0.05
      }
    }
  }
  return(list(loglike=loglike,param=param))
}