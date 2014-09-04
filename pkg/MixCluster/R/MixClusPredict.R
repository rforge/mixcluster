MixClusPredict <- function(x, rule){
  if (!is.matrix(x))
    stop("x must be a matrix")
  
  if (ncol(x)!=rule@data@e)
    stop("x must be composed with the same variables that the variables into rule")
  
  
  data <- MixClusData(x, rep(NA,nrow(x)), rule@data@kind)
  
  model=rule@model$model
  
  data@proba <- MixClusProba(data, rule@param, model)
  model@loglike <- sum(log(rowSums(data@proba)))
  model@bic <- model@loglike - model@nbparam*0.5*log(data@n)
  data@tik <- data@proba/rowSums(data@proba)
  data@partition <- rep(0, data@n)
  for (k in 1:model@g)
    data@partition[which(rowSums(sweep(data@tik, 1, data@tik[,k], "-")<=0)==model@g)] <- k
  
  model@icl <- 0
  for (k in 1:model@g)
    model@icl <- model@icl + sum(log(data@proba[which(data@partition==k),k]))
  model@icl <- model@icl - model@nbparam*0.5*log(data@n)
  
  for (k in 1:model@g){
    data@condexpec[[k]] <- matrix(0, data@n, data@e)
    sup <- matrix(0,data@n,data@e)
    inf <- matrix(0,data@n,data@e)
    for (j in 1:data@e){
      bound <- cbind(rep(-Inf, data@n), rep(Inf, data@n))
      bound[data@o[[j]], ] <- findbounds(data@x[data@o[[j]], j], rule@param@beta[[k]][[j]])
      if (data@kind[j]==1){
        bound[,1] <- bound[,1]-10**(-6) 
      }
      sup[,j] <- bound[,2]
      inf[,j] <- bound[,1]
    }
    
    for (i in 1:data@n)
      data@condexpec[[k]][i,] <- mtmvnorm(mean=rep(0,data@e), sigma=rule@param@correl[[k]], lower=inf[i,], upper=sup[i,])$tmean
  }
  
  return(new("MixClusResults",  model=list(model=model), param=rule@param, data=data, priors=rule@priors, initparam=rule@param))
}