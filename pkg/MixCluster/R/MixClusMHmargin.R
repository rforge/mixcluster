setGeneric ( name= "MixClusMHmargin",  def = function(x, param, prior, latent, correl, j){ standardGeneric("MixClusMHmargin")})


# setMethod( f = "MixClusMHmargin", 
#            signature(x="numeric", param="MixClusParam_continuous", prior="MixClusPrior_continuous", latent="matrix", correl="matrix", j="numeric"), 
#            definition = function(x, param, prior, latent, correl, j){
#              
#              ## sampling of the candidate
#              cand <- rep(0,2)
#              # sampling of the sd
#              ck <- prior@co + 0.5*length(x)
#              if (length(x)>1){
#                Ck <- prior@Co + 0.5*((length(x)-1)*var(x) + length(x)*prior@No/(length(x)+prior@No) * ((mean(x)-prior@bo)**2))
#                
#              }else{
#                Ck <- prior@Co + 0.5*(length(x)*prior@No/(length(x)+prior@No) * ((mean(x)-prior@bo)**2))
#                
#              }
#              
#              
#              cand[2] <- sqrt(1/rgamma(1, ck, Ck)) 
#              # sampling of the mean
#              bk <- prior@bo*prior@No/(length(x)+prior@No) + mean(x)*length(x)/(length(x)+prior@No)
#              Bkcand <- (cand[2]**2)/(length(x)+prior@No)
#              Bk <- (param@sd**2)/(length(x)+prior@No)
#              cand[1] <- rnorm(1, bk, sqrt(Bk))
#              
#              if (ncol(latent)>1)
#                mutilde <-  rowSums(sweep(latent, 2, as.numeric(correl[j, -j] %*% solve(correl[-j, -j]) ),FUN="*"))
#              else
#                mutilde <-  as.numeric(correl[j, -j] %*% solve(correl[-j, -j]))*as.numeric(latent)
#              
#              sigmatilde <- sqrt( 1 - correl[j, -j] %*% solve(correl[-j, -j])  %*%correl[-j, j] )  
#              
#              
#              # ratio of acceptance 
#              rho <- exp( (dgamma((1/cand[2])**2, prior@co, prior@Co, log=TRUE) + dnorm(cand[1], prior@bo, sqrt((cand[2]**2)/prior@No), log=TRUE)
#                           - dgamma((1/cand[2])**2, ck, Ck, log=TRUE) - dnorm(cand[1], bk, sqrt(Bkcand), log=TRUE)
#                           + sum(dnorm((x-cand[1])/cand[2], mutilde, sigmatilde, log=TRUE)) - length(x)*log(cand[2])) -
#                            ( dgamma((1/param@sd)**2, prior@co, prior@Co, log=TRUE) + dnorm(param@mu, prior@bo, sqrt((param@sd**2)/prior@No), log=TRUE)
#                              - dgamma((1/param@sd)**2, ck, Ck, log=TRUE) - dnorm(param@mu, bk, sqrt(Bk), log=TRUE)
#                              + sum(dnorm((x-param@mu)/param@sd, mutilde, sigmatilde, log=TRUE)) - length(x)*log(param@sd)
#                            )
#              )
#              
#              
#              
#              if ((!is.na(rho))&&(runif(1) < rho))
#                param <- new("MixClusParam_continuous", mu=cand[1], sd=cand[2])            
#              
#              return(param)
#            }
# )
# 




setMethod( f = "MixClusMHmargin", 
           signature(x="numeric", param="MixClusParam_continuous", prior="MixClusPrior_continuous", latent="matrix", correl="matrix", j="numeric"), 
           definition = function(x, param, prior, latent, correl, j){
             
             ## sampling of the candidate
             cand <- c(param@mu, param@sd)
             
             if ((runif(1)>0.5)&&(length(x)>1)){
               # sampling of the sd
               ck <- prior@co + 0.5*length(x)
               if (length(x)>1){
                 Ck <- prior@Co + 0.5*((length(x)-1)*var(x) + length(x)*prior@No/(length(x)+prior@No) * ((mean(x)-prior@bo)**2))
                 
               }else{
                 Ck <- prior@Co + 0.5*(length(x)*prior@No/(length(x)+prior@No) * ((mean(x)-prior@bo)**2))
                 
               }
               cand[2] <- sqrt(1/rgamma(1, ck, Ck)) 
               
               
               rho <-  (dgamma((1/cand[2])**2, prior@co, prior@Co, log=TRUE) - dgamma((1/cand[2])**2, ck, Ck, log=TRUE)) -
                             ( dgamma((1/param@sd)**2, prior@co, prior@Co, log=TRUE) - dgamma((1/param@sd)**2, ck, Ck, log=TRUE) )
               
             }else{
               bk <- prior@bo*prior@No/(length(x)+prior@No) + mean(x)*length(x)/(length(x)+prior@No)
               Bk <- (param@sd**2)/(length(x)+prior@No)
               cand[1] <- rnorm(1, bk, sqrt(Bk))
              
               rho <- ( dnorm(cand[1], prior@bo, sqrt((cand[2]**2)/prior@No), log=TRUE) - dnorm(cand[1], bk, sqrt(Bk), log=TRUE))-
                             (dnorm(param@mu, prior@bo, sqrt((param@sd**2)/prior@No), log=TRUE)- dnorm(param@mu, bk, sqrt(Bk), log=TRUE))
                               
             }

             
         
             # sampling of the mean

             
             if (ncol(latent)>1)
               mutilde <-  rowSums(sweep(latent, 2, as.numeric(correl[j, -j] %*% solve(correl[-j, -j]) ),FUN="*"))
             else
               mutilde <-  as.numeric(correl[j, -j] %*% solve(correl[-j, -j]))*as.numeric(latent)
             
             sigmatilde <- sqrt( 1 - correl[j, -j] %*% solve(correl[-j, -j])  %*%correl[-j, j] )  
             
             
             # ratio of acceptance 
             rho <- exp( rho + sum(dnorm((x-cand[1])/cand[2], mutilde, sigmatilde, log=TRUE)) - length(x)*log(cand[2]) -
                           sum(dnorm((x-param@mu)/param@sd, mutilde, sigmatilde, log=TRUE)) - length(x)*log(param@sd))
             
             
             
             if ((!is.na(rho))&&(runif(1) < rho))
               param <- new("MixClusParam_continuous", mu=cand[1], sd=cand[2])            
             
             return(param)
           }
)

setMethod( f = "MixClusMHmargin", 
           signature(x="numeric", param="MixClusParam_integer", prior="MixClusPrior_integer", latent="matrix", correl="matrix", j="numeric"), 
           definition = function(x, param, prior, latent, correl, j){
             
             ## sampling of the candidate
             cand <- new("MixClusParam_integer", beta=rgamma(1, sum(x) + prior@ao ,length(x) + prior@Ao ))         
             
             mutilde <-  rowSums(sweep(latent, 2, as.numeric(correl[j, -j] %*% solve(correl[-j, -j]) ),FUN="*"))
             sigmatilde <- sqrt( 1 - correl[j, -j] %*% solve(correl[-j, -j])  %*%correl[-j, j] )  
             
             bounds <- findbounds(x, param) 
             boundscand <- findbounds(x, cand)                        
             
             # ratio of acceptance 
             rho <-  exp( dgamma(param@beta, sum(x) + prior@ao ,length(x) + prior@Ao, log=TRUE ) - dgamma(cand@beta, sum(x) + prior@ao ,length(x) + prior@Ao, log=TRUE ) 
                          - dgamma(param@beta,  prior@ao , prior@Ao, log=TRUE ) + dgamma(cand@beta, prior@ao , prior@Ao, log=TRUE ) 
                          + sum(log(pnorm(boundscand[,2], mutilde, sigmatilde) - pnorm(boundscand[,1], mutilde, sigmatilde)))  
                          - sum(log(pnorm(bounds[,2], mutilde, sigmatilde) - pnorm(bounds[,1], mutilde, sigmatilde))) )
             
             if ((!is.na(rho))&&(runif(1) < rho))
               param <- cand          
             
             return(param)
           }
)


setMethod( f = "MixClusMHmargin", 
           signature(x="numeric", param="MixClusParam_ordinal", prior="MixClusPrior_ordinal", latent="matrix", correl="matrix", j="numeric"), 
           definition = function(x, param, prior, latent, correl, j){
             
             ## sampling of the candidate      
             val <- prior@a
             for (u in unique(x))
               val[ u+1 ] <- sum(x==u) + prior@a[ u+1 ]
             
             cand <- new("MixClusParam_ordinal", beta=rdiri(val))      
             
             mutilde <-  rowSums(sweep(latent, 2, as.numeric(correl[j, -j] %*% solve(correl[-j, -j]) ),FUN="*"))
             sigmatilde <- sqrt( 1 - correl[j, -j] %*% solve(correl[-j, -j])  %*%correl[-j, j] )  
             
             bounds <- findbounds(x, param) 
             boundscand <- findbounds(x, cand)                        
             
             # ratio of acceptance 
             rho <-  exp( ddiri(param@beta, val, log=TRUE) - ddiri(cand@beta, val, log=TRUE)
                          + ddiri(cand@beta, prior@a, log=TRUE) - ddiri(param@beta, prior@a, log=TRUE)
                          + sum(log(pnorm(boundscand[,2], mutilde, sigmatilde) - pnorm(boundscand[,1], mutilde, sigmatilde))) 
                          -  sum(log(pnorm(bounds[,2], mutilde, sigmatilde) - pnorm(bounds[,1], mutilde, sigmatilde))) )
          
             if ((!is.na(rho))&&(runif(1) < rho))
               param <- cand            
             
             return(param)
           }
)