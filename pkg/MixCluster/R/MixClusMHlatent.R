




# This function performs the sampling of the latent vectors (z and y)
# It uses the first Metropolis-Hastings algorithm described in the article
MixClusMHlatent <- function(data, model, param, latent){
  
  # Sampling of the candidate: zstar is the class membership
  # and ystar is the latent vector related to the copula
  # zstar is uniformly sample
  zstar <- sample(1:model@g, data@n, replace=TRUE)
  
  # initialization of ystar 
  ystar <- matrix(0, data@n, data@e)
  
  # initialization of the vector containing the probability of acceptance
  rho <- rep(1, data@n)
  
  for (j in 1:data@e){
    # bound and bounstar are the bound for y and ystar respectively
    # first column is the lower bound and second column is the upper bound
    bound <- cbind(rep(-Inf, data@n), rep(Inf, data@n))
    boundstar <- bound
    
    for (k in 1:model@g){
      who <- which((!is.na(data@x[,j]))*(zstar==k) ==1)
      if (length(who)>0)
        boundstar[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])      
      
      who <- which((!is.na(data@x[,j]))*(latent@z==k) ==1)
      if (length(who)>0)
        bound[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
      
    }
    # sampling of ystar according to its posterior distribution
    # under the locally independent model
    who <- which(((boundstar[, 1]<Inf)*(boundstar[,2]>-Inf))==1)
    ystar[who, j] <- rtnorm(length(who), lower=boundstar[who, 1], upper=boundstar[who, 2])
    
    # update the acceptance probability according to the proposal part
    if (data@kind[j]!=1)
      rho[who] <- rho[who] * dtnorm(latent@y[who, j], lower=bound[who, 1], upper=bound[who, 2]) / dtnorm(ystar[who, j], lower=boundstar[who, 1], upper=boundstar[who, 2])
    
    if (length(who)<data@n){
      rho[-who] <- 0
    }
  }
  
  # update the acceptance probability according to the stationnary distribution part
  for (k in 1:model@g){
    if (any(zstar==k))
      rho[zstar==k] <- rho[zstar==k] * dmvnorm(ystar[zstar==k,], sigma=param@correl[[k]]) * param@pi[k]
    
    if (any(latent@z==k))
      rho[latent@z==k] <- rho[latent@z==k] / (dmvnorm(latent@y[latent@z==k,], sigma=param@correl[[k]]) * param@pi[k])
    
  }
  
  # simulation to determine which candidates are accepted
  accepte <- which((runif(data@n) < rho)==1)
  
  #update the latent vector
  if (length(accepte)>0){
    latent@z[accepte] <- zstar[accepte]
    latent@y[accepte,] <- ystar[accepte,]
  }
  
  
  return(latent)
  
}



# This function performs the sampling of the latent vectors (z and y)
# It uses the first Metropolis-Hastings algorithm described in the article
MixClusMHlatent_semi_supervized <- function(data, model, param, latent){
  
  # Sampling of the candidate: zstar is the class membership
  # and ystar is the latent vector related to the copula
  # zstar is uniformly sample
  zstar <- sample(1:model@g, data@n, replace=TRUE)
  
  # initialization of ystar 
  ystar <- matrix(0, data@n, data@e)
  
  # initialization of the vector containing the probability of acceptance
  rho <- rep(1, data@n)
  
  for (j in 1:data@e){
    # bound and bounstar are the bound for y and ystar respectively
    # first column is the lower bound and second column is the upper bound
    bound <- cbind(rep(-Inf, data@n), rep(Inf, data@n))
    boundstar <- bound
    
    for (k in 1:model@g){
      who <- which((!is.na(data@x[,j]))*(zstar==k) ==1)
      if (length(who)>0)
        boundstar[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])      
      
      who <- which((!is.na(data@x[,j]))*(latent@z==k) ==1)
      if (length(who)>0)
        bound[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
      
    }
    # sampling of ystar according to its posterior distribution
    # under the locally independent model
    who <- which(((boundstar[, 1]<Inf)*(boundstar[,2]>-Inf))==1)
    ystar[who, j] <- rtnorm(length(who), lower=boundstar[who, 1], upper=boundstar[who, 2])
    
    # update the acceptance probability according to the proposal part
    if (data@kind[j]!=1)
      rho[who] <- rho[who] * dtnorm(latent@y[who, j], lower=bound[who, 1], upper=bound[who, 2]) / dtnorm(ystar[who, j], lower=boundstar[who, 1], upper=boundstar[who, 2])
    
    if (length(who)<data@n){
      rho[-who] <- 0
    }
  }
  
  # update the acceptance probability according to the stationnary distribution part
  for (k in 1:model@g){
    if (any(zstar==k))
      rho[zstar==k] <- rho[zstar==k] * dmvnorm(ystar[zstar==k,], sigma=param@correl[[k]]) * param@pi[k]
    
    if (any(latent@z==k))
      rho[latent@z==k] <- rho[latent@z==k] / (dmvnorm(latent@y[latent@z==k,], sigma=param@correl[[k]]) * param@pi[k])
    
  }
  
  rho[which(!is.na(data@partition))] <- -1
  # simulation to determine which candidates are accepted
  accepte <- which((runif(data@n) < rho)==1)
  
  #update the latent vector
  if (length(accepte)>0){
    latent@z[accepte] <- zstar[accepte]
    latent@y[accepte,] <- ystar[accepte,]
  }
  
  
  return(latent)
  
}