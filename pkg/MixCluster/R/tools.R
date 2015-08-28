# this function samples a vector from a Dirichlet distribution
rdiri <- function(alpha){
  output <- rep(0, length(alpha))
  for (u in 1:length(output))
    output[u] <- rgamma(1, alpha[u], 1)
  
  return(output/sum(output))
}

# this function computes the density of vector according to a Dirichlet distribution
ddiri <- function(x, alpha, log=FALSE){
  output <- sum((alpha-1)*log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
  if (log==FALSE){
    output <- exp(output)
  }
  return(output)  
}


# The function "findbounds" determines the interval of y 
# which provides the same x. Obviously, this function depends about
# the nature of the observed variable x
setGeneric ( name= "findbounds",  def = function(x, param){ standardGeneric("findbounds")})

setMethod( f = "findbounds", 
           signature(x="numeric", param="MixClusParam_continuous"), 
           definition = function(x, param){
             binf <- (x - param@mu)/param@sd
             return(  cbind(binf, binf))
           }
)

setMethod( f = "findbounds", 
           signature(x="numeric", param="MixClusParam_integer"), 
           definition = function(x, param){
             return(cbind(qnorm(ppois(x-1, param@beta)), qnorm(ppois(x, param@beta)))) 
           }
)

setMethod( f = "findbounds", 
           signature(x="numeric", param="MixClusParam_ordinal"), 
           definition = function(x, param){
             val <- c(0, cumsum(param@beta))
             val[ length(param@beta) + 1 ] <- 1
             return(cbind(qnorm(val[ x + 1]), qnorm(val[ x + 2 ])))
           }
)


setGeneric ( name= "Savemargin",  def = function(backup, param, nbiter){ standardGeneric("Savemargin")})

setMethod( f = "Savemargin", 
           signature(backup="MixClusParam_continuous", param="MixClusParam_continuous", nbiter="numeric"), 
           definition = function(backup, param, nbiter){
             backup@mu <- backup@mu + param@mu/nbiter
             backup@sd <- backup@sd + param@sd/nbiter
             return(backup)
           }
)


setMethod( f = "Savemargin", 
           signature(backup="MixClusParam_integer", param="MixClusParam_integer", nbiter="numeric"), 
           definition = function(backup, param, nbiter){
             backup@beta <- backup@beta + param@beta/nbiter
             return(backup)
           }
)


setMethod( f = "Savemargin", 
           signature(backup="MixClusParam_ordinal", param="MixClusParam_ordinal", nbiter="numeric"), 
           definition = function(backup, param, nbiter){
             backup@beta <- backup@beta + param@beta/nbiter
             return(backup)
           }
)


Buildbackup <- function(kind, param){
  backup <- param
  backup@pi <- backup@pi*0
  for (k in 1:length(param@pi)){
    backup@correl[[k]] <- backup@correl[[k]]*0
    for (j in 1:length(backup@beta[[k]])){
      if (kind[j]==1){
        backup@beta[[k]][[j]]@mu <- 0
        backup@beta[[k]][[j]]@sd <- 0
      }else{
        backup@beta[[k]][[j]]@beta <- backup@beta[[k]][[j]]@beta*0
      }
      
    }
  }
  return(backup)
}