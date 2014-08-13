sample_correlations <- function(latent,correlations){
  for (k in unique(latent[,1])){
    who <- which(latent[,1]==k)
    if (length(who)>prod(dim(correlations[[k]]))/2){
      correlations[[k]] <- sample_correlations_class(latent[who,-1])
    }else{
      #correlations[[k]] <- sample_correlations_class(t(as.matrix(latent[who,-1])))
    }
  }
  return(correlations)
}

sample_correlations_homo <- function(latent,correlations){
  tmp <-  sample_correlations_class(latent[,-1])
  for (k in unique(latent[,1])){
    correlations[[k]] <- tmp
  }
  return(correlations)
}


sample_correlations_uncoeff <- function(latent,correlations){
  for (k in unique(latent[,1])){
    who <- which(latent[,1]==k)
    if (length(who)>prod(dim(correlations[[k]]))/2){
      tmp <- 1/rgamma(1,(length(who)+1)/2,length(who)*(sum(cov(latent[,-1])) - sum(diag(cov(latent[who,-1]))))/2)
      tmp <- matrix(tmp,ncol(correlations[[k]]), ncol(correlations[[k]]))
      diag(tmp) <- 1
      if (det(tmp)>0)
        correlations[[k]] <- tmp
    }else{
      #correlations[[k]] <- sample_correlations_class(t(as.matrix(latent[who,-1])))
    }
  }

  
  return(correlations)
}


rW <- function(S0,nu){
  sS0 <- chol(S0)
  Z <- matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
  t(Z)%*%Z
}
sample_correlations_class <- function(y){
  tmp <- solve(rwish(3+0.5*nrow(y),solve(diag(1,ncol(y)) + 0.5 * (nrow(y)-1)* cov(y))))
  #tmp <- solve(rwish(nrow(y),solve( (nrow(y)-1)* cov(y))))
  return(diag(1/sqrt(diag(tmp))) %*% tmp %*% diag(1/sqrt(diag(tmp))))
}
