genere_cand_binaire <- function(x){
  return(rbeta(1,sum(1-x)+1/2,sum(x)+1/2))
}

compute_logproposal_binaire <- function(beta, x, me, sig){
  output <- pnorm(qnorm(beta), me, sig)
  output[which(x==1)] <- 1 - output[which(x==1)]
  output <-  sum(log(output)) + dbeta(beta,1/2,1/2,log=TRUE) - dbeta(beta,sum(1-x)+1/2,sum(x)+1/2,log=TRUE)
  return(output)
}

genere_y_binaire2 <- function(beta, x, me, sig){
  b <- matrix(qnorm(beta),length(x),2)
  b[which(x==0),1] <- -Inf
  b[which(x==1),2] <- Inf
  who <- which(b[,1] < b[,2])
  output <- rep(NA,length(x))
  output[who] <- rtnorm(length(who), me[who], sig, lower=b[who,1], upper=b[who,2])
  return(output)
}

MH_binaire <- function(x,latent,param,j){                     
  for (k in unique(latent[,1])){
    if (sum(latent[,1]==k)>1){
      me <-  rowSums(sweep(as.matrix(latent[which(latent[,1]==k), -c(1,j+1)]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
      sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )
      
      beta_star <- genere_cand_binaire(x[which(latent[,1]==k),j])      
      rho2 <- exp(compute_logproposal_binaire(beta_star, x[which(latent[,1]==k),j], me, sig) - compute_logproposal_binaire(param@margins[[k]][j,1], x[which(latent[,1]==k),j], me, sig))

      if (rho2>runif(1)){
        param@margins[[k]][j,1] <- beta_star
      }
      latent[which(latent[,1]==k),j+1]  <-genere_y_binaire2(param@margins[[k]][j,1], x[which(latent[,1]==k),j], me, sig)
    }

  }
  return(list(latentj=latent[,j+1],param=param))
}