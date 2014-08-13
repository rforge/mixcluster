sample_margins_and_gaussian <- function(x,latent,param,priors,kind){
  for (j in 1:ncol(x)){
    if (kind[j]==1){
      tmp <- MH_continue(x,latent,param,j,priors)
    }else if (kind[j]==2){
      tmp <- MH_entier(x,latent,param,j)
    }else if (kind[j]==3){
      tmp <- MH_binaire(x,latent,param,j)
    }
    latent[,j+1] <- tmp$latentj
    param <- tmp$param
  }
  return(list(latent=latent,param=param))
}
