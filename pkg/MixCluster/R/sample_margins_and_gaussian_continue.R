genere_cand_continue <- function(x,priors){
  cand <- rep(0,2)
  ck <- priors$co + 0.5*length(x)
  Ck <- priors$Co + 0.5*((length(x)-1)*var(x) + length(x)*priors$No/(length(x)+priors$No) * ((mean(x)-priors$bo)**2))
  cand[2] <- sqrt(1/rgamma(1,ck,Ck)) 
  
  bk <- priors$bo*priors$No/(length(x)+priors$No) + mean(x)*length(x)/(length(x)+priors$No)
  Bk <- (cand[2]**2)/(length(x)+priors$No)
  cand[1] <- rnorm(1,bk,sqrt(Bk))
  
  return(cand)
}


compute_logproposal_continue <- function(bet, x, y, me, sig, priors){
  ck <- priors$co + 0.5*length(x)
  Ck <- priors$Co + 0.5*((length(x)-1)*var(x) + length(x)*priors$No/(length(x)+priors$No) * ((mean(x)-priors$bo)**2))
  bk <- priors$bo*priors$No/(length(x)+priors$No) + mean(x)*length(x)/(length(x)+priors$No)
  Bk <- (bet[2]**2)/(length(x)+priors$No)
 # print(c(ck,Ck,bk,Bk))
  output <- dgamma(1/(bet[2]**2),priors$co,priors$Co,log=TRUE) + dnorm(bet[1],priors$bo,sqrt((bet[2]**2)/priors$No), log=TRUE) - dgamma(1/(bet[2]**2),ck,Ck,log=TRUE)- dnorm(bet[1],bk,sqrt(Bk), log=TRUE)
#  print(output)
  output <- output + sum(dnorm(y,me,sig,log=TRUE)) - length(x)*log(bet[2])
#  print(output)
  return(output)
}


MH_continue <- function(x,latent,param,j,priors){
  priors <- priors[[j]]
  for (k in unique(latent[,1])){
    if (sum(latent[,1]==k)>2){
      me <-  rowSums(sweep(as.matrix(latent[which(latent[,1]==k), -c(1,j+1)]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
      sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )
      
      beta_star <- genere_cand_continue(x[which(latent[,1]==k),j],priors)
      rho2 <- exp( compute_logproposal_continue(beta_star, x[which(latent[,1]==k),j], genere_y_continue(x[which(latent[,1]==k),j],beta_star), me, sig, priors) - compute_logproposal_continue(param@margins[[k]][j,], x[which(latent[,1]==k),j], genere_y_continue(x[which(latent[,1]==k),j],param@margins[[k]][j,]), me, sig, priors) )
      
        if (rho2>runif(1)){
        param@margins[[k]][j,] <- beta_star
      }
      latent[which(latent[,1]==k),j+1] <- genere_y_continue(x[which(latent[,1]==k),j],param@margins[[k]][j,])
    }
  }
  return(list(latentj=latent[,j+1],param=param))
}