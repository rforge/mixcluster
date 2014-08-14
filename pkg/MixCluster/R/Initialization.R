##########################################################################

## Initialization of the parameters

GCMMparam_init_EM_loc_indpt <- function(g,x,kind){
  pi <- as.numeric(rdiri(rep(2,g)))
  margins <- list()
  gamma <- list()
  for (k in 1:g){
    if (all(kind!=4)){
      margins[[k]] <- matrix(0,ncol(x),2)
    }else{
      margins[[k]] <- matrix(0,ncol(x),max(x[,which(kind==4)])+1)
    }
    for (j in 1:ncol(x)){
      if (kind[j]==1){
        margins[[k]][j,1:2] <- c(x[sample(1:nrow(x),1),j],sd(x[,j]))
      }else if (kind[j]==2){
        margins[[k]][j,1:2] <- c(x[sample(1:nrow(x),1),j]+0.2,0)
      }else if (kind[j]>2){
        margins[[k]][j,1:length(unique(x[,j]))] <- (rdiri(rep(1/2, length(unique(x[,j])))))
      }
    }
    gamma[[k]] <- diag(1,ncol(x),ncol(x))
  }
  return(GCMMparam(pi,margins,gamma))
}

GCMMparam_init <- function(g,x,kind){
  base <- XEM_loc_indpt(x,g,kind)
  for (u in 1:15){
    cand <- XEM_loc_indpt(x,g,kind)
    if (cand$loglike>base$loglike)
      base <- cand
  }
  return(base$param)
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
    }else if (kind[j]>2){
      for (k in 1:g){
        for (u in unique(x[,j])){
          who <- which(x[,j]==u)
          proba[who,k] <- proba[who,k] * param@margins[[k]][j,u+1]
        }

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
       }else if (kind[j]>2){
          for (u in unique(x[,j])){
            who <- which(x[,j]==u)
            param@margins[[k]][j,u+1] <- sum(proba[who,k])/sum(proba[,k])
            new_proba[who,k] <- new_proba[who,k] * param@margins[[k]][j,u+1]
          }
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
  return(list(loglike=loglike,param=param))
}


##########################################################################

## Initialization of the latent variables

GCMMlatent_init <- function(x,g,param,kind){
  latent <- matrix(1,nrow(x),ncol(x)+1)
  if (g>1){
    proba <- matrix(param@proportions,nrow(x),g,byrow=TRUE)
    for (j in 1:ncol(x)){
      for (k in 1:g)
        proba[,k] <- proba[,k] * pdfmargin(x[,j],param@margins[[j]][[k]])
    }
    latent[,1] <- rowSums( sweep( t( apply( as.matrix( proba/rowSums(proba) ),1,cumsum ) ), 1, runif( nrow(x)) , "<=") )+1
  }
  for (k in unique(latent[,1])){
    for (j in 1:ncol(x))
      latent[which(latent[,1]==k),j+1] <- sample_yj_indpt(x[which(latent[,1]==k),j],param@margins[[j]][[k]])
  }
  return(latent)
}
