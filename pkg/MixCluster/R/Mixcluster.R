MixClusanalysis <- function(x, g, kind=NULL, model="free", nbchains=1, burn=100, iterGibbs=1000, param=NULL, latent=NULL){
  cat("  Initialization of the Gibbs sampler\n")
  
  if (missing(x))
    stop("Data set is missing!")
  
  if (any(is.na(x)))
    stop("Data set containing missing values cannot be analyzed")
  
  if (missing(g))
    stop("Class number is missing!")
  
  # Define the nature of the variables if it is missing.
  if (is.null(kind)){
    kind <- rep(1,ncol(x))
    for (j in 1:ncol(x)){
      if (all(x[,j]==ceiling(x[,j]))){
        if (max(x[,j])==1){
          kind[j] <- 3
        }else if (length(unique(x[,j]))<5){
          kind[j] <- 4
        }else{
          kind[j] <- 2
        }
      }
    }
  }else if ( (length(kind)!=ncol(x)) || (is.numeric(kind)==FALSE)){
      stop("The nature of the variables is badly specified")
  }
  
  if (all(kind==3))
    cat("  The identifiability of the model is not ensured\n")
  
  if (is.null(param)){
    param <- GCMMparam_init(g,x,kind)
  }else if (class(param)!="MixClusparam"){
    cat(" warning: the parameters wanted to be used for the initialization were not taken into account since they are not a instance of the MixClusparam class!")
    param <- GCMMparam_init(g,x,kind)
  }
  
  tmp_param <- param
  param@margins <- list()
  for (j in 1:ncol(x)){
    param@margins[[j]] <- list()
    if (kind[j]==1){for (k in 1:g){ param@margins[[j]][[k]] <- GCMMparam_continue( tmp_param@margins[[k]][j, 1:2] ) }}
    if (kind[j]==2){for (k in 1:g){ param@margins[[j]][[k]] <- GCMMparam_entier( tmp_param@margins[[k]][j, 1] ) }}
    if (kind[j]==3){for (k in 1:g){ param@margins[[j]][[k]] <- GCMMparam_ordinal( tmp_param@margins[[k]][j, 1:length(unique(x[,j]))] )}}
  }
  
  rm(tmp_param)
  
  if (is.null(latent)){
    latent <- GCMMlatent_init(x, g, param, kind)
  }else if ((is.matrix(latent)==FALSE) || (nrow(latent)!=nrow(x)) || ((ncol(latent)-1) != ncol(x))){
    cat(" warning: the latent variables used to initialize the initialization were not taken into account since 
        they are not on the good format: a matrix where the first column denotes the class membership of the individuals and 
        where the other columns denote the latent Gaussian vector.")
    latent <- GCMMlatent_init(x,g,param,kind)
  }
  
   
  
  if (model!="all"){
    output <- GCMM(x, g, kind, burn, iterGibbs, param, latent, model)
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- GCMM(x, g, kind, burn, iterGibbs, param, latent, model)
        if (cand@criteria$loglike>output@criteria$loglike)
          output <- cand
      }
    }
  }else if (model=="all"){
    output <- list(free=GCMM(x, g, kind, burn, iterGibbs, param, latent, "free"), homo=GCMM(x, g, kind, burn, iterGibbs, param, latent, "homo"), indpt=GCMM(x, g, kind, burn, iterGibbs, param, latent, "indpt"))
    if (nbchains>1){
      for (u in 2:nbchains){
        cand <- list(free=GCMM(x, g, kind, burn, iterGibbs, param, latent, "free"), homo=GCMM(x, g, kind, burn, iterGibbs, param, latent, "homo"), indpt=GCMM(x, g, kind, burn, iterGibbs, param, latent, "indpt"))
        if (cand$free@criteria$loglike>output$free@criteria$loglike)
          output$free <- cand$free
        if (cand$homo@criteria$loglike>output$homo@criteria$loglike)
          output$homo <- cand$homo
        if (cand$indpt@criteria$loglike>output$indpt@criteria$loglike)
          output$indpt <- cand$indpt
      }
    }
  }else{
    stop(" Model is not specified")
  }
  return(output)
}



GCMM<-function(x, g, kind, burn, iterGibbs, param, latent, model){
  if (model=="free"){
    model <- GCMM_model_free(model)
  }else if (model=="homo"){
    model <- GCMM_model_homo(model)
  }else if (model=="indpt"){
    model <- GCMM_model_indpt(model)
  }
  # initialization of the objects of GCM
   priors <- GCMMpriors(x, kind)
   cat("  Burn-in of the free model\n")
   for (it in 1:burn){
     # sampling of the latent variables 
     latent <- MH1_latent(x, g, latent, param, kind)
    # sampling of the marings and gaussian parameters
     for (k in unique(latent[,1])){
       if (sum(latent[,1]==k)>2){
         for (j in 1:ncol(x)){
           me <-  rowSums(sweep(as.matrix(latent[which(latent[,1]==k), -c(1,j+1)]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
           sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )  
           beta_star <- candidate_maring_sampler(x[which(latent[,1]==k),j], param@margins[[j]][[k]], priors[[j]])
           rho2 <- exp( logproba_param_MH2(x[which(latent[,1]==k),j], beta_star, me, sig, priors[[j]])-logproba_param_MH2(x[which(latent[,1]==k),j], param@margins[[j]][[k]], me, sig, priors[[j]]))    
 
           
           if (rho2>runif(1)){
             param@margins[[j]][[k]] <- beta_star
             if (class(param@margins[[j]][[k]])=="MixClusparam_continue")
               latent[which(latent[,1]==k),j+1] <- (x[which(latent[,1]==k),j] - param@margins[[j]][[k]]@param[1] )/ param@margins[[j]][[k]]@param[2]
           }
           if (class(param@margins[[j]][[k]])!="MixClusparam_continue")
             latent[which(latent[,1]==k),j+1] <- genere_y_cond(x[which(latent[,1]==k),j],param@margins[[j]][[k]], me, sig)
           
         }
       }
     }
    
    # sampling of the correlation matrices
    param@correlations <- sample_correlations(latent,param@correlations,model)
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
  }
  
  save_param <- param
  
  cat("  Parameter inference of the free model\n")
  for (it in 1:iterGibbs){
    # sampling of the latent variables 
    latent <- MH1_latent(x, g, latent, param, kind)
    
    # sampling of the marings and gaussian parameters
    for (k in unique(latent[,1])){
      if (sum(latent[,1]==k)>2){
        for (j in 1:ncol(x)){
          me <-  rowSums(sweep(as.matrix(latent[which(latent[,1]==k), -c(1,j+1)]), 2, as.numeric(param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j]) ),FUN="*"))
          sig <- sqrt( param@correlations[[k]][j,j]- param@correlations[[k]][j,-j] %*% solve(param@correlations[[k]][-j,-j])  %*%param@correlations[[k]][-j,j] )  
          beta_star <- candidate_maring_sampler(x[which(latent[,1]==k),j], param@margins[[j]][[k]], priors[[j]])
          rho2 <- exp( logproba_param_MH2(x[which(latent[,1]==k),j], beta_star, me, sig, priors[[j]])-logproba_param_MH2(x[which(latent[,1]==k),j], param@margins[[j]][[k]], me, sig, priors[[j]]))    
          if (rho2>runif(1)){
            param@margins[[j]][[k]] <- beta_star
            if (class(param@margins[[j]][[k]])=="MixClusparam_continue")
              latent[which(latent[,1]==k),j+1] <- (x[which(latent[,1]==k),j] - param@margins[[j]][[k]]@param[1] )/ param@margins[[j]][[k]]@param[2]
          }
          if (class(param@margins[[j]][[k]])!="MixClusparam_continue")
            latent[which(latent[,1]==k),j+1] <- genere_y_cond(x[which(latent[,1]==k),j],param@margins[[j]][[k]], me, sig)
          
          save_param@margins[[j]][[k]]@param <- save_param@margins[[j]][[k]]@param + param@margins[[j]][[k]]@param
        }
      }
    }
    
    # sampling of the correlation matrices
    param@correlations <- sample_correlations(latent,param@correlations,model)
    for (k in 1:g)
      save_param@correlations[[k]] <- save_param@correlations[[k]] + param@correlations[[k]]
    
    # sampling of the vector of proportions
    param@proportions <- sample_proportions(table(c(1:g,latent[,1])))
    save_param@proportions <- save_param@proportions + param@proportions
  }
  
   

   
  cat("  Results bulding\n")
  save_param <- normsave(save_param,iterGibbs)
   passe <- save_param
   passe@margins<- list()
   for (k in 1:g){
     passe@margins[[k]] <- list()
     if (any(kind>2)){
       passe@margins[[k]] <- matrix(0,ncol(x),max(x[which(kind>2),])+1)
     }else{
       passe@margins[[k]] <- matrix(0,ncol(x),2)
     }
     for (j in 1:ncol(x)){
       passe@margins[[k]][j,1:length(save_param@margins[[j]][[k]]@param)] <- save_param@margins[[j]][[k]]@param
     }
   }
   
   criteria <- compute_criteria(x,passe,kind, model)
   
   output <- MixClusresultsctr(param=passe, tik=criteria$tik, partition=criteria$z,
                               criteria=list(loglike=criteria$loglike, bic=criteria$bic, icl=criteria$icl),
                               data=list(nrow=nrow(latent), ncol=length(kind), kind=kind, data=x),
                               current_latent=latent,model=model@name, expect_y=list(), test_expect_y=rep(0,g))
   nam <- colnames(x)
   if (is.null(nam))
     nam <- paste("X",1:ncol(x),sep="")
   

   for (k in 1:g){
     colnames(output@param@margins[[k]])=paste("param",1:ncol(output@param@margins[[k]]),sep="")
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
