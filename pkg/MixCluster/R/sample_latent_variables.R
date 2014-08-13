# candidat generation
genere_candidat_latent <- function(x,g,param,kind){
  latent_cand <- cbind(sample(1:g,nrow(x),replace=TRUE),matrix(0,nrow(x),ncol(x)))
  for (k in unique(latent_cand[,1])){
    for (j in 1:ncol(x)){
      latent_cand[which(latent_cand[,1]==k),j+1] <- sample_yj_indpt(x[which(latent_cand[,1]==k),j],kind[j],param@margins[[k]][j,])
    }
  }
  return(latent_cand)
}

compute_logproposal_latent <- function(x,g,latent,param,kind){
  proposal <- log(param@proportions[latent[,1]])
  for (k in 1:g){
    if (any(latent[,1]==k)){
      proposal[which(latent[,1]==k)] <- proposal[which(latent[,1]==k)] + dmvnorm(latent[which(latent[,1]==k),-1],mean=rep(0,ncol(latent)-1),sigma=param@correlations[[k]],log=TRUE)
      for (j in which(kind>1)){
        if (kind[j]==2){
          who <- which((latent[,1]==k)*(is.na(latent[,j+1])==FALSE) == 1)    
          binf <- qnorm(ppois(x[who,j]-1,param@margins[[k]][j,1]))
          bsup <- qnorm(ppois(x[who,j],param@margins[[k]][j,1]))
          
          proposal[who] <- proposal[who]  - dtnorm(latent[who,j+1],lower=binf,upper=bsup,log=TRUE)
          if (any((latent[,1]==k)*(is.na(latent[,j+1])) == 1)){
            who <- which((latent[,1]==k)*(is.na(latent[,j+1])) == 1)
            proposal[who] <- -Inf
          }
          
        }else if (kind[j]==3){
          who <- which((latent[,1]==k)*(is.na(latent[,j+1])==FALSE)*(x[,j]==0) == 1)    
          proposal[who] <- proposal[who]  - dtnorm(latent[who,j+1],lower=-Inf,upper=qnorm(param@margins[[k]][j,1]),log=TRUE)
          who <- which((latent[,1]==k)*(is.na(latent[,j+1])==FALSE)*(x[,j]==1) == 1)    
          proposal[who] <- proposal[who]  - dtnorm(latent[who,j+1],lower=qnorm(param@margins[[k]][j,1]),upper=Inf,log=TRUE)
          
          if (any((latent[,1]==k)*(is.na(latent[,j+1])) == 1)){
            who <- which((latent[,1]==k)*(is.na(latent[,j+1])) == 1)
            proposal[who] <- -Inf
          }
        }
      }
    }
  }
  return(proposal)
}

sample_latent_variables <- function(x,g,latent,param,kind){
  l_c <- genere_candidat_latent(x,g,param,kind)
  pro <- exp(compute_logproposal_latent(x,g,l_c,param,kind) - compute_logproposal_latent(x,g,latent,param,kind))
  pro[is.na(pro)] <- 0
  who <- which(runif(nrow(latent))<pro)
  latent[who,] <- l_c[who,]
  return(latent)
}


genere_y_entier <- function(x,margins){
  output <- rep(0,length(x))
  for (u in unique(x)){
    if ( (qnorm(ppois(u-1,margins[1]))<Inf) && (qnorm(ppois(u,margins[1]))>-Inf) ){
      output[which(x==u)] <- rtnorm(sum(x==u),0,1,lower=qnorm(ppois(u-1,margins[1])),upper=qnorm(ppois(u,margins[1])))
    }else if (qnorm(ppois(u-1,margins[1]))==Inf){
      output[which(x==u)] <- NA
    }else{
      output[which(x==u)] <- NA
      
    }
  }
  return(output)
}

logproba_y_entier <- function(x,y,margins){
  res <- 0
  for (u in unique(x)){
    if (qnorm(ppois(u-1,margins[1]))<qnorm(ppois(u,margins[1]))){
      res <- res + sum(dtnorm(y[which(x==u)],0,1,lower=qnorm(ppois(u-1,margins[1])),upper=qnorm(ppois(u,margins[1])),log=TRUE))
    }else{ 
      res <- -Inf
      break
    }
  }
  return(res)
}

genere_y_binaire <- function(x,margins){
  output <- rep(0,length(x))
  if (any(x==0)){
    output[which(x==0)] <- rtnorm(sum(x==0),0,1,lower=-Inf,upper=qnorm(margins[1]))
  }  
  if (any(x==1)){
    output[which(x==1)] <- rtnorm(sum(x==1),0,1,lower=qnorm(margins[1]),upper=Inf)
  }
  return(output)
}

logproba_y_binaire <- function(x,y,margins){
  sum(dtnorm(y[which(x==0)],0,1,lower=-Inf,upper=qnorm(margins),log=TRUE)) + sum(dtnorm(y[which(x==1)],0,1,lower=qnorm(margins),upper=Inf,log=TRUE))
}

genere_y_continue <- function(x,margins){
  (x-margins[1])/margins[2]
}

sample_yj_indpt<-function(x,kind,margins){
  if (kind==1){
    output <- genere_y_continue(x,margins)
  }else if (kind==2){
    output <- genere_y_entier(x,margins)
  }else if (kind==3){
    output <- genere_y_binaire(x,margins)
  }
  return(output)
}