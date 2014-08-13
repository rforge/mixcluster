genere_cand_entier <- function(x){
  return(rgamma(1,sum(x)+1/2,length(x)+1/2))
}


compute_logproposal_entier <- function(beta, x, me, sig){
  output <- sum(log(pnorm(qnorm(ppois(x,beta)), me, sig) - pnorm(qnorm(ppois(x-1,beta)), me, sig))) + dgamma(beta,1/2,1/2,log=TRUE) - dgamma(beta,sum(x)+1/2,length(x)+1/2,log=TRUE)
  return(output)
}


genere_y_entier2 <- function(beta, x, me, sig){
  b <- cbind(qnorm(ppois(x-1, beta)), qnorm(ppois(x, beta)))
  who <- which(b[,1] < b[,2])
  output <- rep(NA,length(x))
  output[who] <- rtnorm(length(who), me[who], sig, lower=b[who,1], upper=b[who,2])
  return(output)
}

MH_entier <- function(x,latent,param,j){                      
  for (k in unique(latent[,1])){
    if (sum(latent[,1]==k)>1){
      me <-  rowSums(sweep(as.matrix(latent[which(latent[,1]==k), -c(1,j+1)]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
      sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )
      
      beta_star <- genere_cand_entier(x[which(latent[,1]==k),j])
      rho2 <- exp(compute_logproposal_entier(beta_star, x[which(latent[,1]==k),j], me, sig ) - compute_logproposal_entier(param@margins[[k]][j,1], x[which(latent[,1]==k),j], me, sig))
                  
      if ((is.na(rho2)==FALSE)&&(rho2>runif(1))){
        param@margins[[k]][j,1] <- beta_star
      }
  
      latent[which(latent[,1]==k),j+1] <- genere_y_entier2(param@margins[[k]][j,1], x[which(latent[,1]==k),j], me, sig)
    }
  }
  return(list(latentj=latent[,j+1],param=param))
}