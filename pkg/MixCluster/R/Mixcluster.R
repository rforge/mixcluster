MixClusanalysis <- function(x, g, kind=NULL, model="free", nbchains=1, burn=100, iterGibbs=1000, param=NULL, latent=NULL){
  cat("  Initialization of the Gibbs sampler\n")
  
  # Define the nature of the variables if it is missing.
  if (is.null(kind)){
    kind <- rep(1,ncol(x))
    for (j in 1:ncol(x)){
      if (all(x[,j]==ceiling(x[,j]))){
        if (max(x[,j])==1){
          kind[j] <- 3
        }else{
          kind[j] <- 2
        }
      }
    }
  }
  if (all(kind==3))
    cat("  The identifiability of the model is not ensured\n")
  
  if (is.null(param))
    param <- GCMMparam_init(g,x,kind)
  
  if (is.null(latent))
    latent <- GCMMlatent_init(x,g,param,kind)
  
  if (model=="free"){
    output <- GCMM_free(x, g, kind, burn, iterGibbs, param, latent)
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- GCMM_free(x, g, kind, burn, iterGibbs, param, latent)
        if (cand@criteria$loglike>output@criteria$loglike)
          output <- cand
      }
    }
  }else if (model=="homo"){
    output <- GCMM_homo(x, g, kind, burn, iterGibbs, param, latent)
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- GCMM_homo(x, g, kind, burn, iterGibbs, param, latent)
        if (cand@criteria$loglike>output@criteria$loglike)
          output <- cand
      }
    }
  }else if (model=="indpt"){
    output <- GCMM_indpt(x, g, kind, burn, iterGibbs, param, latent)
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- GCMM_indpt(x, g, kind, burn, iterGibbs, param, latent)
        if (cand@criteria$loglike>output@criteria$loglike)
          output <- cand
      }
    }
  }else if (model=="all"){
    output <- list(free=GCMM_free(x, g, kind, burn, iterGibbs, param, latent), homo=GCMM_homo(x, g, kind, burn, iterGibbs, param, latent), indpt=GCMM_indpt(x, g, kind, burn, iterGibbs, param, latent))
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- list(free=GCMM_free(x, g, kind, burn, iterGibbs, param, latent), homo=GCMM_homo(x, g, kind, burn, iterGibbs, param, latent), indpt=GCMM_indpt(x, g, kind, burn, iterGibbs, param, latent))
        if (cand$free@criteria$loglike>output$free@criteria$loglike)
          output$free <- cand$free
        if (cand$homo@criteria$loglike>output$homo@criteria$loglike)
          output$homo <- cand$homo
        if (cand$indpt@criteria$loglike>output$indpt@criteria$loglike)
          output$indpt <- cand$indpt
      }
    }
  }else{
    cat("  Error: model is not specified\n")
  }
  return(output)
}



GCMM_free<-function(x, g, kind, burn, iterGibbs, param, latent){
  # initialization of the objects of GCM
   priors <- GCMMpriors(x,kind)
   cat("  Burn-in of the free model\n")
   for (it in 1:burn){
     # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param    
    # sampling of the correlation matrices
    param@correlations <- sample_correlations(latent,param@correlations)
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
  }
  
  save_param <- param
  
  cat("  Parameter inference of the free model\n")
  for (it in 1:iterGibbs){
    # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param    
    # sampling of the correlation matrices
    param@correlations <- sample_correlations(latent,param@correlations)
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
    
    save_param <- addparam(save_param,param)
  }
  
  cat("  Results bulding\n")
  save_param <- normsave(save_param,iterGibbs)
  criteria <- compute_criteria(x,save_param,kind)
  
   
   output <- MixClusresultsctr(param=save_param, tik=criteria$tik, partition=criteria$z,
                               criteria=list(loglike=criteria$loglike, bic=criteria$bic, icl=criteria$icl),
                               data=list(nrow=nrow(latent), ncol=length(kind), kind=kind, data=x),
                               current_latent=latent,model="free", expect_y=list(), test_expect_y=rep(0,g))
   nam <- colnames(x)
   if (is.null(nam))
     nam <- paste("X",1:ncol(x),sep="")
   

   for (k in 1:g){
     colnames(output@param@margins[[k]])=paste("param",1:g,sep="")
     rownames(output@param@margins[[k]])=nam
     colnames(output@param@correlations[[k]])=nam
     rownames(output@param@correlations[[k]])=nam
   }
   names(output@param@correlations)=paste("Class",1:g,sep="")
   names(output@param@margins)=paste("Class",1:g,sep="")
   colnames(output@tik)=paste("Class",1:g,sep="")
   cat("\n")
  return(output)
}



GCMM_homo<-function(x, g, kind, burn, iterGibbs, param, latent){
  priors <- GCMMpriors(x,kind)
  cat("  Burn-in of the homoscedastic model\n")
  for (it in 1:burn){
    # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param
    # sampling of the correlation matrices
    param@correlations <- sample_correlations_homo(latent,param@correlations) 
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
  }
  
  save_param <- param
  
  cat("  Parameter inference of the homoscedastic model\n")
  for (it in 1:iterGibbs){
    # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param
    # sampling of the correlation matrices
    param@correlations <- sample_correlations_homo(latent,param@correlations) 
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
    
    save_param <- addparam(save_param,param)
  }
  
  cat("  Results building\n")
  save_param <- normsave(save_param,iterGibbs)
  
  criteria <- compute_criteria(x,save_param,kind,hetero=FALSE)
  
  output <- MixClusresultsctr(param=save_param, tik=criteria$tik, partition=criteria$z,
                              criteria=list(loglike=criteria$loglike, bic=criteria$bic, icl=criteria$icl),
                              data=list(nrow=nrow(latent), ncol=length(kind), kind=kind, data=x),
                              current_latent=latent,
                              model="homo", expect_y=list(), test_expect_y=rep(0,g))
  nam <- colnames(x)
  if (is.null(nam))
    nam <- paste("X",1:ncol(x),sep="")
  
  
  for (k in 1:g){
    colnames(output@param@margins[[k]])=paste("param",1:g,sep="")
    rownames(output@param@margins[[k]])=nam
    colnames(output@param@correlations[[k]])=nam
    rownames(output@param@correlations[[k]])=nam
  }
  names(output@param@correlations)=paste("Class",1:g,sep="")
  names(output@param@margins)=paste("Class",1:g,sep="")
  colnames(output@tik)=paste("Class",1:g,sep="")
  cat("\n")
  return(output)
}



GCMM_indpt<-function(x, g, kind, burn, iterGibbs, param, latent){
  priors <- GCMMpriors(x,kind)
  
  cat("  Burn-in of the locally independent model\n")
  
  for (it in 1:burn){
    # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
  }
  
  save_param <- param
  save_latent <- latent
  
  cat("  Parameter inference of the locally independent model\n")
  for (it in 1:iterGibbs){
    # sampling of the latent variables 
    latent <- sample_latent_variables(x,g,latent,param,kind)
    # sampling of the marings and gaussian parameters
    tmp <- sample_margins_and_gaussian(x,latent,param,priors,kind)
    latent <- tmp$latent
    param <- tmp$param
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))

    save_param <- addparam(save_param,param)
  }
  
  cat("  Results building\n")
  save_param <- normsave(save_param,iterGibbs)
  
  criteria <- compute_criteria(x,save_param,kind,hetero=FALSE,cim=TRUE)
  
  
  output <- MixClusresultsctr(param=save_param, tik=criteria$tik, partition=criteria$z,
                              criteria=list(loglike=criteria$loglike, bic=criteria$bic, icl=criteria$icl),
                              data=list(nrow=nrow(latent), ncol=length(kind), kind=kind, data=x),
                              current_latent=latent,model="indpt", expect_y=list(), test_expect_y=rep(0,g))
  nam <- colnames(x)
  if (is.null(nam))
    nam <- paste("X",1:ncol(x),sep="")
  
  
  for (k in 1:g){
    colnames(output@param@margins[[k]])=paste("param",1:g,sep="")
    rownames(output@param@margins[[k]])=nam
    colnames(output@param@correlations[[k]])=nam
    rownames(output@param@correlations[[k]])=nam
  }
  names(output@param@correlations)=paste("Class",1:g,sep="")
  names(output@param@margins)=paste("Class",1:g,sep="")
  colnames(output@tik)=paste("Class",1:g,sep="")
  cat("\n")
  return(output)
}
