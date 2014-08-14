# candidat generation
setGeneric ( name= "logproba_q_yj",  def = function(x, y, param){ standardGeneric("logproba_q_yj")}) 
setMethod( f = "logproba_q_yj", signature(x="numeric", y="numeric", param="MixClusparam_entier"), definition = function(x, y, param){dnorm(y,log=TRUE) - log(ppois(x,param@param)-ppois(x-1,param@param))})
setMethod( f = "logproba_q_yj", signature(x="numeric", y="numeric", param="MixClusparam_ordinal"), definition = function(x, y, param){dnorm(y,log=TRUE) - log(param@param[x+1]) })


logproba_param_MH1 <- function(x, g, latent, param, kind){
  proposal <- log(param@proportions[latent[,1]])
  for (k in 1:g){
    if (any(latent[,1]==k)){
      proposal[which(latent[,1]==k)] <- proposal[which(latent[,1]==k)] + dmvnorm(latent[which(latent[,1]==k),-1],mean=rep(0,ncol(latent)-1),sigma=param@correlations[[k]],log=TRUE)
      for (j in which(kind>1)){
        who <- which((latent[,1]==k)*(is.na(latent[,j+1])==FALSE) == 1)
        proposal[who] <- proposal[who] - logproba_q_yj(x[who,j], latent[who,j+1], param@margins[[j]][[k]])     
        
        if (any((latent[,1]==k)*(is.na(latent[,j+1])) == 1)){
         who <- which((latent[,1]==k)*(is.na(latent[,j+1])) == 1)
         proposal[who] <- -Inf
        }
        
      }
    }
  }
  return(proposal)
}

MH1_latent <- function(x, g, latent, param, kind){
  latent_cand <- cbind(sample(1:g,nrow(x),replace=TRUE),matrix(0,nrow(x),ncol(x)))
  for (j in 1:ncol(x)){
    for (k in unique(latent_cand[,1]))
      latent_cand[which(latent_cand[,1]==k),j+1] <- sample_yj_indpt(x[which(latent_cand[,1]==k),j], param@margins[[j]][[k]])
  }
  
  pro <- exp(logproba_param_MH1(x,g,latent_cand ,param,kind) - logproba_param_MH1(x,g,latent,param,kind))
  pro[is.na(pro)] <- 0
  
  who <- which(runif(nrow(latent))<pro)
  latent[who,] <- latent_cand[who,]
  
  return(latent)
}




