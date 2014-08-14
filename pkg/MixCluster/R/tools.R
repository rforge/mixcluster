rdiri <- function(alpha){
  output <- rep(0,length(alpha))
  for (u in 1:length(output))
    output[u] <- rgamma(1,alpha[u],1)
  
  return(output/sum(output))
}

ddiri <- function(x, alpha, loga=FALSE){
  output <- sum((alpha-1)*log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
  if (loga==FALSE){
    output <- exp(output)
  }
  return(output)  
}

setClass(Class="GCMM_model_free", representation=representation(name="character"), prototype = prototype(name=character(0)))
setClass(Class="GCMM_model_homo", representation=representation(name="character"), prototype = prototype(name=character(0)))
setClass(Class="GCMM_model_indpt", representation=representation(name="character"), prototype = prototype(name=character(0)))
GCMM_model_free <- function(param){ new("GCMM_model_free",name=param)}
GCMM_model_homo  <- function(param){ new("GCMM_model_homo",name=param)}
GCMM_model_indpt <- function(param){ new("GCMM_model_indpt",name=param)}

setClass(Class="MixClusparam_continue", representation=representation(param="numeric"), prototype = prototype(param=numeric(0)))
setClass(Class="MixClusparam_entier"  , representation=representation(param="numeric"), prototype = prototype(param=numeric(0)))
setClass(Class="MixClusparam_ordinal" , representation=representation(param="numeric"), prototype = prototype(param=numeric(0)))
GCMMparam_continue<- function(param){ new("MixClusparam_continue",param=param)}
GCMMparam_entier  <- function(param){ new("MixClusparam_entier",param=param)}
GCMMparam_ordinal <- function(param){ new("MixClusparam_ordinal",param=param)}

setGeneric ( name= "pdfmargin",  def = function(x, param){ standardGeneric("pdfmargin")}) 
setMethod( f = "pdfmargin", signature(x="numeric", param="MixClusparam_continue"), definition = function(x, param){dnorm(x, param@param[1], param@param[2])} )
setMethod( f = "pdfmargin", signature(x="numeric", param="MixClusparam_entier"), definition = function(x, param){dpois(x, param@param)} )
setMethod( f = "pdfmargin", signature(x="numeric", param="MixClusparam_ordinal"), definition = function(x, param){param@param[x+1] })



setGeneric ( name= "sample_yj_indpt",  def = function(x, param){ standardGeneric("sample_yj_indpt")}) 
setMethod( f = "sample_yj_indpt", signature(x="numeric", param="MixClusparam_continue"), definition = function(x, param){(x-param@param[1])/param@param[2] })
  
setMethod( f = "sample_yj_indpt", signature(x="numeric", param="MixClusparam_entier"),
           definition = function(x, param){
            output <- rep(0,length(x))
            for (u in unique(x)){
              if ( (qnorm(ppois(u-1,param@param))<Inf) && (qnorm(ppois(u,param@param))>-Inf) ){
                output[which(x==u)] <- rtnorm(sum(x==u),0,1,lower=qnorm(ppois(u-1,param@param)),upper=qnorm(ppois(u,param@param)))
              }else if (qnorm(ppois(u-1,param@param))==Inf){
                output[which(x==u)] <- NA
              }else{
                output[which(x==u)] <- NA
                
              }
            }
            return(output)
          }
)


setMethod( f = "sample_yj_indpt", signature(x="numeric", param="MixClusparam_ordinal"), definition = function(x, param){
  output <- rep(0,length(x))
  tmpval <- cumsum(c(0,param@param))
  for (u in unique(x)){
    output[which(x==u)] <- rtnorm(sum(x==u), 0, 1, lower=tmpval[u+1], upper=tmpval[u+2])
  }
  return(output)
}
)



GCMMpriors <- function(x,kind){
  priors <- list();
  for (j in 1:ncol(x)){
    priors[[j]]=list()
  
    if (kind[j]){
      priors[[j]] <- list()
      priors[[j]]$bo <- mean(x[,j])
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
          }else if (kind[j]==4){
            binf[cp] <- qnorm(cumsum(c(0,param@margins[[k]][j,]))[(x[i,j]+1)])
            bsup[cp] <- qnorm(cumsum(param@margins[[k]][j,])[(x[i,j]+1)])
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

compute_criteria <- function(x,param,kind, model){
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
  if (model@name=="free"){
      nbparam <- (g-1)+g*(d*(d-1)/2 ) + 2*g*sum(kind==1) + g*sum(kind==2)
      if (any(kind>2)){
        for (j in which(kind>2))
        nbparam <- nbparam + g*max(x[,j])
      }
  }else if (model@name=="homo"){
        nbparam <- (g-1)+ d*(d-1)/2 + 2*g*sum(kind==1) + g*sum(kind==2)
        if (any(kind>2)){
          for (j in which(kind>2))
            nbparam <- nbparam + g*max(x[,j])
        }
  }else if (model@name=="indpt"){
      nbparam <- (g-1) + 2*g*sum(kind==1) + g*sum(kind==2)
      if (any(kind>2)){
        for (j in which(kind>2))
          nbparam <- nbparam + g*max(x[,j])
      }
  }else{
    print("pb nb param")
  }
  print(nbparam)
  bic <- loglike - nbparam*0.5*log(nrow(x))
  icl <- icl - nbparam*0.5*log(nrow(x))
  return(list(tik=tik, loglike=loglike, bic=bic, icl=icl, nbparam=nbparam, z=z))
}


addparam <- function(save_param,param){
  save_param@proportions <- save_param@proportions + param@proportions
  for (k in 1:length(save_param@proportions)){
    for (j in 1:length(save_param@margins))
      save_param@margins[[j]][[k]] <- save_param@margins[[j]][[k]] + param@margins[[j]][[k]]
    save_param@correlations[[k]] <- save_param@correlations[[k]] + param@correlations[[k]]
  }
  return(save_param)
}

normsave <- function(param,itermax){
  param@proportions <- param@proportions/(itermax+1)
  for (k in 1:length(param@proportions)){
    for (j in 1:length(param@margins))
    param@margins[[j]][[k]]@param <- param@margins[[j]][[k]]@param/(itermax+1)
    param@correlations[[k]] <- param@correlations[[k]]/(itermax+1)
  }
  return(param)
}