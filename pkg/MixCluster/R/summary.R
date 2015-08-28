setMethod(
  f="summary",
  signature=c("MixClusResults"),
  function(object, ...){
      cat("**************************************************************\n\n")
      cat("DATA SET description:","\n")
      cat("* Number of individuals: ", object@data@n,"\n")
      cat("* Number of variables (continuous, integer, ordinal): ", paste(sum(object@data@kind==1), 
                                                                          sum(object@data@kind==2),
                                                                          sum(object@data@kind==3),
                                                                          sep=", "),"\n\n")
      cat("**************************************************************\n\n")
      cat("MODEL description:","\n")
      cat("* Number of classes:   ", length(object@param@pi),"\n")
      cat("* Model:            ", class(object@model$model),"\n")
      cat("* Log-likelihood:  ", object@model$model@loglike,"\n")
      cat("* BIC:             ", object@model$model@bic,"\n")
      cat("* ICL:             ", object@model$model@icl,"\n")
      
      cat("**************************************************************\n\n")
      cat("PARAMETERS description:","\n")
      for (k in 1:length(object@param@pi)){
        cat("* Class", k,"\n")         
        cat("* Proportion            ", object@param@pi[k],"\n")
        if (any(object@data@kind==1)){
          cat("* Margins parameters of the continuous variables (Gaussian distribution):\n")
          tmp <- matrix(0, sum(object@data@kind==1), 2)
          rownames(tmp)=colnames(object@data@x)[which(object@data@kind==1)]
          colnames(tmp) <- c("Mean","Sd")
          cp <- 1
          for (j in which(object@data@kind==1)){
            tmp[cp, ] <- c(object@param@beta[[k]][[j]]@mu, object@param@beta[[k]][[j]]@sd)  
            cp <- cp+1
          }

          print(tmp, digits=2)
          cat("\n")
        }
        if (any(object@data@kind==2)){
          cat("* Margins parameters of the integer variables (Poisson distribution):\n")
   
          
          tmp <- matrix(0, sum(object@data@kind==2), 1)
          rownames(tmp)=colnames(object@data@x)[which(object@data@kind==2)]
          colnames(tmp) <- c("Mean")
          cp <- 1
          for (j in which(object@data@kind==2)){
            tmp[cp, ] <- c(object@param@beta[[k]][[j]]@beta)  
            cp <- cp+1
          }
          print(tmp, digits=2)
          cat("\n")
        }
        if (any(object@data@kind==3)){
          cat("* Margins parameters of the ordinal variables (multinomial distribution):\n")
          lg <- 0
          for (j in which(object@data@kind==3))
            lg <- max(lg, max(object@data@x[object@data@o[[j]],j])+1)
          
          tmp <- matrix(0, sum(object@data@kind==3), lg)
          rownames(tmp) <- colnames(object@data@x)[which(object@data@kind==3)]
          colnames(tmp) <- paste("Prob.",0:(lg-1),sep="")
           cp <- 1
          for (j in which(object@data@kind==3)){
            tmp[cp, 1:(max(object@data@x[object@data@o[[j]],j])+1)] <- object@param@beta[[k]][[j]]@beta
            cp <- cp + 1
          }
          print(tmp, digits=2)
          cat("\n")
        }
        cat("\n")
        if (class(object@model$model)=="MixClusModel_hetero"){
          cat("* Correlation matrix:\n")
          print(object@param@correl[[k]],digits=2)
        }
        cat("\n")
        
        cat("*****************************\n")
          
      }
      if (class(object@model$model)!="MixClusModel_hetero"){
        cat("* Correlations matrix:\n")
        print(object@param@correl[[k]],digits=2)
        
      }

    return(invisible())
  }
)