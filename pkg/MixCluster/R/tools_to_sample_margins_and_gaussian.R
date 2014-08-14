setGeneric ( name= "candidate_maring_sampler",  def = function(x, param, priors){ standardGeneric("candidate_maring_sampler")}) 


setMethod( f = "candidate_maring_sampler", signature(x="numeric", param="MixClusparam_continue", priors="list"), definition = function(x, param, priors){
    cand <- rep(0,2)
    ck <- priors$co + 0.5*length(x)
    Ck <- priors$Co + 0.5*((length(x)-1)*var(x) + length(x)*priors$No/(length(x)+priors$No) * ((mean(x)-priors$bo)**2))
    cand[2] <- sqrt(1/rgamma(1,ck,Ck)) 
    
    bk <- priors$bo*priors$No/(length(x)+priors$No) + mean(x)*length(x)/(length(x)+priors$No)
    Bk <- (cand[2]**2)/(length(x)+priors$No)
    cand[1] <- rnorm(1,bk,sqrt(Bk))
    return(GCMMparam_continue(cand))
  }
)

setMethod( f = "candidate_maring_sampler", signature(x="numeric", param="MixClusparam_entier", priors="list"), definition = function(x, param, priors){
    GCMMparam_entier(rgamma(1,sum(x)+1/2,length(x)+1/2))
  }
)

setMethod( f = "candidate_maring_sampler", signature(x="numeric", param="MixClusparam_ordinal", priors="list"), definition = function(x, param, priors){
    val <- rep(0,length(param@param))
    for (u in unique(x))
      val[u+1] <- sum(x==u)+1/2
    
    GCMMparam_ordinal(rdiri(val))
  }
)


setGeneric ( name= "logproba_param_MH2",  def = function(x, param, me, sig, priors){ standardGeneric("logproba_param_MH2")}) 


setMethod( f = "logproba_param_MH2", signature(x="numeric", param="MixClusparam_continue", me="numeric", sig="matrix",  priors="list"), 
           definition = function(x,  param, me, sig, priors){
             ck <- priors$co + 0.5*length(x)
             Ck <- priors$Co + 0.5*((length(x)-1)*var(x) + length(x)*priors$No/(length(x)+priors$No) * ((mean(x)-priors$bo)**2))
             bk <- priors$bo*priors$No/(length(x)+priors$No) + mean(x)*length(x)/(length(x)+priors$No)
             Bk <- (param@param[2]**2)/(length(x)+priors$No)
             
             output <- dgamma(1/(param@param[2]**2),priors$co,priors$Co,log=TRUE) + dnorm(param@param[1],priors$bo,sqrt((param@param[2]**2)/priors$No), log=TRUE) -
               dgamma(1/(param@param[2]**2),ck,Ck,log=TRUE)- dnorm(param@param[1],bk,sqrt(Bk), log=TRUE) + sum(dnorm((x-param@param[1])/param@param[2],me,sig,log=TRUE)) - length(x)*log(param@param[2])
             
             return(output)
           }
)



setMethod( f = "logproba_param_MH2", signature(x="numeric", param="MixClusparam_entier", me="numeric", sig="matrix",  priors="list"), 
           definition = function(x,  param, me, sig, priors){
             sum(log(pnorm(qnorm(ppois(x,param@param)), me, sig) - pnorm(qnorm(ppois(x-1,param@param)), me, sig))) + dgamma(param@param,1/2,1/2,log=TRUE) - dgamma(param@param,sum(x)+1/2,length(x)+1/2,log=TRUE)
           }
)




setMethod( f = "logproba_param_MH2", signature(x="numeric", param="MixClusparam_ordinal", me="numeric", sig="matrix",  priors="list"), 
           definition = function(x,  param, me, sig, priors){
             
             val <- rep(0,length(param@param))
             for (u in unique(x))
               val[u+1] <- sum(x==u)+1/2
             
             tmpval <- c(-Inf,qnorm(cumsum(param@param)[-length(param@param)]), Inf)
             
             output <-  log(ddiri(param@param,rep(1/2,length(param@param)))) - log(ddiri(param@param, val)) +             
                         sum(log( pnorm(tmpval[x+2] , me, sig) - pnorm(tmpval[x+1] , me, sig)))
             
             return(output)
           }
)
           

setGeneric ( name= "genere_y_cond",  def = function(x, param, me, sig){ standardGeneric("genere_y_cond")}) 


setMethod( f = "genere_y_cond", signature(x="numeric", param="MixClusparam_entier", me="numeric", sig="matrix"), 
           definition = function(x,  param, me, sig){
             b <- cbind(qnorm(ppois(x-1, param@param)), qnorm(ppois(x, param@param)))
             who <- which((b[,2] - b[,1] )>10^(-4))
             output <- rep(10^8,length(x))
             output[who] <- rtnorm(length(who), me[who], sig, lower=b[who,1], upper=b[who,2])
             return(output)
           }
)



setMethod( f = "genere_y_cond", signature(x="numeric", param="MixClusparam_ordinal", me="numeric", sig="matrix"), 
           definition = function(x,  param, me, sig){
             output <- rep(10^8,length(x))
             tmpval <- c(-Inf,qnorm(cumsum(param@param)[-length(param@param)]), Inf)
             b <- cbind(tmpval[x+1], tmpval[x+2])
             who <- which((b[,2] - b[,1] )>10^(-4))
             output[who] <- rtnorm(length(who),me[who], sig, lower=b[who,1], upper=b[who,2])

             return(output)
           }
)



    
    
