
rW <- function(S0,nu){
  sS0 <- chol(S0)
  Z <- matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
  t(Z)%*%Z
}
sample_correlations_class <- function(y){
  tmp <- solve(rwish(3+0.5*nrow(y),solve(diag(1,ncol(y)) + 0.5 * (nrow(y)-1)* cov(y))))
  #tmp <- solve(rwish(nrow(y),solve( (nrow(y)-1)* cov(y))))
  return(diag(1/sqrt(diag(tmp))) %*% tmp %*% diag(1/sqrt(diag(tmp))))
}

setGeneric ( name= "sample_correlations",  def = function(latent, correlations, model){ standardGeneric("sample_correlations")}) 


setMethod( f = "sample_correlations", signature(latent="matrix", correlations="list", model="GCMM_model_free"),
           definition = function(latent, correlations, model){
             for (k in unique(latent[,1])){
               who <- which(latent[,1]==k)
               if (length(who)>prod(dim(correlations[[k]]))/2){
                 correlations[[k]] <- sample_correlations_class(latent[who,-1])
               }
             }
             return(correlations)
           }
)


setMethod( f = "sample_correlations", signature(latent="matrix", correlations="list", model="GCMM_model_homo"),
           definition = function(latent, correlations, model){
             tmp <-  sample_correlations_class(latent[,-1])
             for (k in unique(latent[,1])){
               correlations[[k]] <- tmp
             }
             return(correlations)
           }
)


setMethod( f = "sample_correlations", signature(latent="matrix", correlations="list", model="GCMM_model_indpt"),
           definition = function(latent, correlations, model){
             return(correlations)
           }
)




