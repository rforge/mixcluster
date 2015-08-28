# Constructor of the class MixClusData
# if kind is not specified, the function automatically detect the 
# nature of each variable.
MixClusData <- function(x, partition, kind=NULL){
  if (is.matrix(x)==FALSE)
    stop("The data set must be a matrix")
  
  for (j in 1:ncol(x)){
    if (is.numeric(x)==FALSE)
      stop("The columns of the data set must be numeric")
  }
  
  if (is.null(kind)){
    kind <- rep(1,ncol(x))
    for (j in 1:ncol(x)){
      who <- which(!is.na(x[,j]))
      if (all(x[who, j]==ceiling(x[who, j]))){
        if (length(unique(x[who,j]))<=5){
          kind[j] <- 3
          x[who,j] <- as.numeric(as.factor(x[who,j])) - 1
        }else{
          kind[j] <- 2
        }
      }
    }
  }else if ( (length(kind)!=ncol(x)) || (is.numeric(kind)==FALSE)){
    stop("The nature of the variables is badly specified")
  }
  o <- list()
  for (j in 1:ncol(x))
    o[[j]] <- which(!is.na(x[,j]))
  
  return(new("MixClusData", x=x, o=o, e=ncol(x), n=nrow(x), kind=kind, partition=as.numeric(partition)))
}


MixClusParam_random <- function(g, data){
  pi <- rdiri(alpha=rep(1/2, g))
  beta <- list()
  correl <- list()
  for (k in 1:g){
    beta[[k]] <- list()
    for (j in 1:data@e){
      if (data@kind[j]==1)
        beta[[k]][[j]] <- new("MixClusParam_continuous", mu=sample(min(data@x[data@o[[j]],j]) : max(data@x[data@o[[j]],j]), 1) ,sd=sd(data@x[data@o[[j]],j])*1.5)
      if (data@kind[j]==2)
        beta[[k]][[j]] <- new("MixClusParam_integer", beta= 1+ sample(min(data@x[data@o[[j]],j]) : max(data@x[data@o[[j]],j]),1))
      if (data@kind[j]==3)
        beta[[k]][[j]] <- new("MixClusParam_ordinal", beta=rdiri(rep(1/2, length(unique(data@x[data@o[[j]], j])))))
    }
    correl[[k]] <- diag(1,data@e,data@e)
  }
  new("MixClusParam", pi=pi, beta=beta, correl=correl)
}

setGeneric ( name= "MixClusProbaMargin",  def = function(x, param){ standardGeneric("MixClusProbaMargin")})

setMethod( f = "MixClusProbaMargin", 
           signature(x="numeric", param="MixClusParam_continuous"), 
           definition = function(x, param)
             return(dnorm(x, param@mu, param@sd))
)


setMethod( f = "MixClusProbaMargin", 
           signature(x="numeric", param="MixClusParam_integer"), 
           definition = function(x, param)
             return(dpois(x, param@beta))
)


setMethod( f = "MixClusProbaMargin", 
           signature(x="numeric", param="MixClusParam_ordinal"), 
           definition = function(x, param)
             return(param@beta[x+1])
)



setGeneric ( name= "MixClusUpdateMargin",  def = function(x, tik, param){ standardGeneric("MixClusUpdateMargin")})

setMethod( f = "MixClusUpdateMargin", 
           signature(x="numeric", tik="numeric", param="MixClusParam_continuous"), 
           definition = function(x, tik, param){
             param@mu <- sum(x*tik)/sum(tik)
             param@sd <- sqrt(sum(((x-param@mu)**2)*tik)/sum(tik))
             return(param)
           }
)


setMethod( f = "MixClusUpdateMargin", 
           signature(x="numeric", tik="numeric", param="MixClusParam_integer"), 
           definition = function(x, tik, param){
             param@beta=sum(x*tik)/sum(tik)
             return(param)
           }
)


setMethod( f = "MixClusUpdateMargin", 
           signature(x="numeric", tik="numeric", param="MixClusParam_ordinal"), 
           definition = function(x, tik, param){
             for (u in unique(x))
               param@beta[u+1] <- sum(tik[which(x==u)])/sum(tik)
             return(param)
           }
)



XEMLocIndpt <- function(data, g){
  loglike <- NA
  while(is.na(loglike)){
    param <- MixClusParam_random(g, data)
    proba <- matrix(param@pi, data@n, length(param@pi), byrow=TRUE)
    for (j in 1:data@e){
      for (k in 1:ncol(proba)){
        proba[data@o[[j]],k] <- proba[data@o[[j]],k] * MixClusProbaMargin(data@x[data@o[[j]], j], param@beta[[k]][[j]])
      }
    }
    loglike <- sum(log(rowSums(proba)))
  }
  
  

  prec <- -Inf
  
  while (loglike-prec>0.01){
    tik <- proba/rowSums(proba)
    param@pi <- colSums(tik)/data@n
    proba <- matrix(param@pi, data@n, length(param@pi), byrow=TRUE)
    for (j in 1:data@e){
      for (k in 1:ncol(proba)){
        param@beta[[k]][[j]] <- MixClusUpdateMargin(data@x[data@o[[j]],j], tik[data@o[[j]],k], param@beta[[k]][[j]])
        proba[data@o[[j]],k] <- proba[data@o[[j]],k] * MixClusProbaMargin(data@x[data@o[[j]], j], param@beta[[k]][[j]])
      }
    }
    prec <- loglike
    loglike <- sum(log(rowSums(proba)))
    if (is.na(loglike))
      loglike <- -Inf
  }
  
  return(list(param=param, loglike=loglike))
}


MixClusParam <- function(data, g, param=NULL){
  if (is.null(param)){
    if (any(is.na(data@partition))){
      param <- XEMLocIndpt(data, g)
      for (u in 2:15){
        cand <- XEMLocIndpt(data, g)
        if (cand$loglike>param$loglike)
          param <- cand
      }
      param <- param$param
    }else{ #if (all(!is.na(data@partition))){
      param <- MixClusParam_random(g, data)
      for (k in 1:g){
        param@pi <- sum(data@partition == k) / data@n
        for (j in 1:data@e)
          param@beta[[k]][[j]] <- MixClusUpdateMargin(data@x[data@o[[j]],j], as.numeric((data@partition==k)[data@o[[j]]]), param@beta[[k]][[j]])
      }
    }
  }
  return(param)
}




setGeneric ( name= "Findprior",  def = function(x, param){ standardGeneric("Findprior")})

setMethod( f = "Findprior", 
           signature(x="numeric",  param="MixClusParam_continuous"), 
           definition = function(x,  param)
             
             new("MixClusPrior_continuous", co=1.28, Co=0.36*var(x), bo=mean(x), No=2.6/(max(x)-min(x)))
)


setMethod( f = "Findprior", 
           signature(x="numeric",  param="MixClusParam_integer"), 
           definition = function(x,  param)
             new("MixClusPrior_integer", ao=1, Ao=1/mean(x))
)


setMethod( f = "Findprior", 
           signature(x="numeric",  param="MixClusParam_ordinal"), 
           definition = function(x, param)
             new("MixClusPrior_ordinal", a=rep(1/2,length(unique(x))))
)


MixClusLatent <- function(data, param){
  if (length(param@pi)>1){
    proba <- matrix(param@pi, data@n, length(param@pi), byrow=TRUE)
    for (j in 1:data@e){
      for (k in 1:ncol(proba)){
        proba[data@o[[j]],k] <- proba[data@o[[j]],k] * MixClusProbaMargin(data@x[data@o[[j]], j], param@beta[[k]][[j]])
      }
    }
    proba <- proba/rowSums(proba)
    proba  <- t(apply(proba, 1, cumsum)) 
    tmp <- matrix(runif(data@n),data@n, ncol(proba))
    z <- ncol(proba) + 1 - rowSums(tmp < proba)
  }else{
    z <- rep(1,data@n)
  }
  
  if (any(!is.na(data@partition)))
    z[which(!is.na(data@partition))] <- data@partition[which(!is.na(data@partition))] 
  
  y <- matrix(0, data@n, data@e)

  for (j in 1:data@e){
    
    # bound and bounstar are the bound for y and ystar respectively
    # first column is the lower bound and second column is the upper bound
    bound <- cbind(rep(-Inf, data@n), rep(Inf, data@n))
    
    for (k in 1:length(param@pi)){
      who <- which((!is.na(data@x[,j]))*(z==k)  ==1)
      if (length(who)>0)
        bound[who, ] <- findbounds(data@x[who, j], param@beta[[k]][[j]])
      
    }
    # sampling of ystar according to its posterior distribution
    # under the locally independent model
    y[, j] <- rtnorm(data@n, lower=bound[, 1], upper=bound[, 2])
  }
  return(new("MixClusLatent", z=z, y=y))
}