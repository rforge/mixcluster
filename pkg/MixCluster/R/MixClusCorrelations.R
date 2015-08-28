
sample_correlations_class <- function(y){
  tmp <- solve( rwish(3+0.5*nrow(y),solve(diag(1,ncol(y)) + 0.5 * (nrow(y)-1)* cov(y))))
  return(diag(1/sqrt(diag(tmp))) %*% tmp %*% diag(1/sqrt(diag(tmp))))
}

setGeneric ( name= "sample_correlations",  def = function(latent, correlations, model){ standardGeneric("sample_correlations")}) 


setMethod( f = "sample_correlations", signature(latent="MixClusLatent", correlations="list", model="MixClusModel_hetero"),
           definition = function(latent, correlations, model){
             for (k in unique(latent@z)){
               if (sum(latent@z==k)>prod(dim(correlations[[k]]))/2){
                 correlations[[k]] <- sample_correlations_class(latent@y[latent@z==k,])
               }
             }
             return(correlations)
           }
)


setMethod( f = "sample_correlations", signature(latent="MixClusLatent", correlations="list", model="MixClusModel_homo"),
           definition = function(latent, correlations, model){
             tmp <-  sample_correlations_class(latent@y)
             for (k in 1:length(correlations)){
               correlations[[k]] <- tmp
             }
             return(correlations)
           }
)


setMethod( f = "sample_correlations", signature(latent="MixClusLatent", correlations="list", model="MixClusModel_indpt"),
           definition = function(latent, correlations, model){
             return(correlations)
           }
)



