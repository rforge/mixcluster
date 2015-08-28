
MixClusGibbs_fix_partition <- function(data, model, param, latent, priors, nbiter, backup){
  
  for (it in 1:nbiter){    
    # Sampling of the margin parameters and sampling of the continuous latent vector
    # The margin parameters are sampled via a Metropolis-Hastings algorithm
    # The continuous vector is sampled from its posterior distribution conditionally on the 
    # margin parameter.
    #print("margin")
    for (j in 1:data@e){
      bounds <-  cbind(rep(-Inf, data@n), rep(Inf, data@n)) 
      for (k in 1:model@g){  
        who <- which((!is.na(data@x[,j]))*(latent@z==k) ==1)
        if (length(who)>1){  
          param@beta[[k]][[j]] <- MixClusMHmargin(data@x[who , j], param@beta[[k]][[j]], priors[[j]], as.matrix(latent@y[who, -j]), param@correl[[k]], j)
          bounds[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
        }else if (length(who)==1){
          #param@beta[[k]][[j]] <- MixClusMHmargin(data@x[who , j], param@beta[[k]][[j]], priors[[j]], matrix(latent@y[who, -j], nrow=1, ncol=data@n-1), param@correl[[k]], j)
          bounds[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
        }
        who <- which(latent@z==k)
        if (length(who)>1)
          mutilde <-  rowSums(sweep(as.matrix(latent@y[who, -j]), 2, as.numeric(param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j]) ),FUN="*"))
        else
          mutilde <- latent@y[who, -j] * (as.numeric(param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j]) ))
        
        sigmatilde <- sqrt( 1 - param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j])  %*%param@correl[[k]][-j, j] )  
        
        latent@y[who, j] <- rtnorm(length(who), mean=mutilde, sd=sigmatilde, lower=bounds[who,1], upper=bounds[who,2])
      }  
    }
    # Sampling of the proportions
    param@pi <- rdiri(table(c(1:model@g,latent@z))-0.5)
    
    # Sampling of the correlation matrices
    param@correl <- sample_correlations(latent, param@correl, model)
    
    # Backup
    if (!missing(backup)){
      backup@pi <- backup@pi + param@pi/nbiter
      for (k in 1:model@g){
        backup@correl[[k]] <- backup@correl[[k]] + param@correl[[k]]/nbiter
        for (j in 1:data@e)
          backup@beta[[k]][[j]] <- Savemargin(backup@beta[[k]][[j]], param@beta[[k]][[j]], nbiter)
      }
    }
  }
  
  if (!missing(backup))
    output <- list(param=param, backup=backup)
  else
    output<- list(param=param)
  
  return(output)
}



One_Gibbs_fix_partition <- function(data, model, param, priors, latent, burn_in, nbiter){
  initparam <- param
  # Burn-in
  if (burn_in>0)
    param <- MixClusGibbs_fix_partition(data, model, param, latent, priors, burn_in)$param
  
  # iteration
  if (nbiter>0)
    param <- MixClusGibbs_fix_partition(data, model, param, latent, priors,  nbiter, Buildbackup(data@kind, param))$backup
  
  # Outputs
  data@partition <- latent@z
  data@proba <- MixClusProba(data, param, model)
  model@loglike <- 0
  data@tik <- matrix(0, data@n, model@g)
  for (k in unique(data@partition)){
    model@loglike <- model@loglike + sum(log(data@proba[which(latent@z==k), k]))
    data@tik[which(latent@z==k),k] <- 1
  }
      
  model@bic <- model@loglike - model@nbparam*0.5*log(data@n)
  model@icl <- model@bic
  
  return(new("MixClusResults",  model=list(model=model), param=param, data=data, priors=priors, initparam=initparam))
}


MixClusGibbs_semi_supervized <- function(data, model, param, latent, priors, nbiter, backup){
  
  for (it in 1:nbiter){
    # Sampling of both latent vectors (z and y)
    # by using a Metropolis-Hastings algorithm
    #print("latent")
    latent <- MixClusMHlatent_semi_supervized(data, model, param, latent)
    
    
    # Sampling of the margin parameters and sampling of the continuous latent vector
    # The margin parameters are sampled via a Metropolis-Hastings algorithm
    # The continuous vector is sampled from its posterior distribution conditionally on the 
    # margin parameter.
    #print("margin")
    for (j in 1:data@e){
      bounds <-  cbind(rep(-Inf, data@n), rep(Inf, data@n)) 
      for (k in 1:model@g){  
        who <- which(((!is.na(data@x[,j]))*(latent@z==k) )==1)
        if (length(who)>1){  
          param@beta[[k]][[j]] <- MixClusMHmargin(data@x[who , j], param@beta[[k]][[j]], priors[[j]], as.matrix(latent@y[who, -j]), param@correl[[k]], j)
          bounds[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
        }else if (length(who)==1){
          #param@beta[[k]][[j]] <- MixClusMHmargin(data@x[who , j], param@beta[[k]][[j]], priors[[j]], matrix(latent@y[who, -j], nrow=1, ncol=data@n-1), param@correl[[k]], j)
          bounds[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
        }
        who <- which(latent@z==k)
        if (length(who)>1)
          mutilde <-  rowSums(sweep(as.matrix(latent@y[who, -j]), 2, as.numeric(param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j]) ),FUN="*"))
        else
          mutilde <- latent@y[who, -j] * (as.numeric(param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j]) ))
        
        sigmatilde <- sqrt( 1 - param@correl[[k]][j, -j] %*% solve(param@correl[[k]][-j, -j])  %*%param@correl[[k]][-j, j] )  
        latent@y[who, j] <- rtnorm(length(who), mean=mutilde, sd=sigmatilde, lower=bounds[who,1], upper=bounds[who,2])
      }  
    }
    
    # Sampling of the proportions
    param@pi <- rdiri(table(c(1:model@g,latent@z)) - 0.5)
    
    # Sampling of the correlation matrices
    param@correl <- sample_correlations(latent, param@correl, model)
    
    # Backup
    if (!missing(backup)){
      backup@pi <- backup@pi + param@pi/nbiter
      for (k in 1:model@g){
        backup@correl[[k]] <- backup@correl[[k]] + param@correl[[k]]/nbiter
        for (j in 1:data@e)
          backup@beta[[k]][[j]] <- Savemargin(backup@beta[[k]][[j]], param@beta[[k]][[j]], nbiter)
      }
      
    }
    
  }
  if (!missing(backup))
    output <- list(param=param, backup=backup)
  else
    output<- list(param=param)
  return(output)
}

One_Gibbs_semi_supervized <- function(data, model, param, priors, latent, burn_in, nbiter){
  initparam <- param
  # Burn-in
  if (burn_in>0)
    param <- MixClusGibbs_semi_supervized(data, model, param, latent, priors, burn_in)$param
  
  # iteration
  if (nbiter>0)
    param <- MixClusGibbs_semi_supervized(data, model, param, latent, priors,  nbiter, Buildbackup(data@kind, param))$backup
  
  # Outputs
  data@proba <- MixClusProba(data, param, model)
  known <- which(!is.na(data@partition))
  unknown <- which(is.na(data@partition))
  
  model@loglike <- 0
  data@tik <- matrix(0, data@n, model@g)
  for (k in unique(data@partition[known])){
    model@loglike <- model@loglike + sum(data@proba[which(data@partition==k), k])
    data@tik[which(data@partition==k), k] <- 1
  }
  
  data@tik[unknown,] <- (data@proba/rowSums(data@proba))[unknown,]
  model@loglike <- model@loglike + sum(log(rowSums(data@proba))[unknown])
  
  model@bic <- model@loglike - model@nbparam*0.5*log(data@n)
  

  data@partition <- rep(0, data@n)
  for (k in 1:model@g)
    data@partition[which(rowSums(sweep(data@tik, 1, data@tik[,k], "-")<=0)==model@g)] <- k
  
  model@icl <- model@bic
  for (k in 1:model@g){
    model@icl <- model@icl + sum(log(data@proba[which(data@partition==k),k]))
  }    
  model@icl <- model@icl - model@nbparam*0.5*log(data@n)
  
  return(new("MixClusResults",  model=list(model=model), param=param, data=data, priors=priors, initparam=initparam))
}


MixClusClassification <- function(x, z, model="hetero", kind=NULL, nbalgo=3, burn_in=300, nbiter=300,  param){
  # Initialization
  g <- length(unique(z[which(!is.na(z))]))
 
  data <- MixClusData(x, z, kind)
  if (missing(param)){
    param <- MixClusParam(data, g)
  }

  
  priors <- MixClusPrior(data, param)
  
  for (k in 1:g)
    data@condexpec[[k]] <- matrix(0,0,0)
  
  model <- MixClusModel(model, g, data, "classification")

  
  if ((g>1)&&(any(is.na(z)))){
    output <- One_Gibbs_semi_supervized(data, model, param, priors, MixClusLatent(data, param), burn_in, nbiter)
    if (nbalgo>1){
      for (it in 2:nbalgo){
        cand <-  One_Gibbs_semi_supervized(data, model, param, priors, MixClusLatent(data, param), burn_in, nbiter)
        if (cand@model$model@bic> output@model$model@bic)
          output <- cand
      }
    }
  }else{
    output <- One_Gibbs_fix_partition(data, model, param, priors, MixClusLatent(data, param), burn_in, nbiter)
    if (nbalgo>1){
      for (it in 2:nbalgo){
        cand <-  One_Gibbs_fix_partition(data, model, param, priors, MixClusLatent(data, param), burn_in, nbiter)
        if (cand@model$model@bic> output@model$model@bic)
          output <- cand
      }
    }
  }
  
  for (k in 1:model@g){
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
      inf[,j] <- bound[,1]
    }
    
    for (i in 1:output@data@n)
      output@data@condexpec[[k]][i,] <- mtmvnorm(mean=rep(0,output@data@e), sigma=output@param@correl[[k]], lower=inf[i,], upper=sup[i,])$tmean
  }
  
  
  return(output)
}
