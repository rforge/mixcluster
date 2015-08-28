MixClusProba <- function(data, param, model){
  proba <- matrix(param@pi, data@n, model@g, byrow=TRUE)
  for (k in 1:model@g){  
    who <- which(data@kind==1)
    
    mu <- rep(0, length(who))
    std <- rep(0, length(who))
    cp <- 1
    for (j in who){
      mu[cp] <- param@beta[[k]][[j]]@mu
      std[cp] <- param@beta[[k]][[j]]@sd
      cp <- cp + 1
    }
    
    sigm <- as.matrix(param@correl[[k]][who,who])
    for (i in 1:data@n){
      coordind <- which((!is.na(data@x[i,who]))==1)
      if (length(coordind)>0){
        proba[i,k] <- proba[i,k] * dmvnorm((data@x[i, which(data@kind==1)][coordind] - mu[coordind])/std[coordind], sigma=as.matrix(sigm[coordind, coordind]) )/ prod(std[coordind]) 
        if (any(data@kind>1)){
          if (any(!is.na(data@x[i,which(data@kind>1)]))){
            vcont <- which(((!is.na(data@x[i,]))*(data@kind==1))==1)
            vdisc <- which(((!is.na(data@x[i,]))*(data@kind>1))==1)
            me <-  as.numeric( t( (param@correl[[k]][vdisc , vcont] %*% solve(param@correl[[k]][vcont, vcont]) ) %*% ((data@x[i, which(data@kind==1)][coordind] - mu[coordind])/std[coordind]) ))
            sig <- as.matrix( param@correl[[k]][vdisc , vdisc] - (param@correl[[k]][vdisc , vcont] %*% solve(param@correl[[k]][vcont, vcont]) ) %*% param@correl[[k]][vcont , vdisc] )
            bounds <- rep(0, 2*length(vdisc))
            tmp1 <- 1
            tmp2 <- tmp1 + length(vdisc)
            for (j in vdisc){
              bounds[c(tmp1,tmp2)] <- findbounds(data@x[i,j], param@beta[[k]][[j]])
              tmp1 <- tmp1+1
              tmp2 <- tmp2+1
            }
            proba[i,k] <- proba[i,k] * pmvnorm(bounds[1:length(vdisc)], bounds[-c(1:length(vdisc))], mean=me, sigma=sig) 
          }
        }
      }else{
        if (any(data@kind>1)){
          if (any(!is.na(data@x[i,which(data@kind>1)]))){
            vdisc <- which(((!is.na(data@x[i,]))*(data@kind>1))==1)
            me <-  rep(0, length(vdisc))
            sig <- as.matrix(param@correl[[k]][vdisc , vdisc])
            bounds <- rep(0, 2*length(vdisc))
            tmp1 <- 1
            tmp2 <- tmp1 + length(vdisc)
            for (j in vdisc){
              bounds[c(tmp1,tmp2)] <- findbounds(data@x[i,j], param@beta[[k]][[j]])
              tmp1 <- tmp1+1
              tmp2 <- tmp2+1
            }
            proba[i,k] <- proba[i,k] * pmvnorm(bounds[1:length(vdisc)], bounds[-c(1:length(vdisc))], mean=me, sigma=sig) 
          }
        }        
      } 
    }
  }
  return(proba)
}