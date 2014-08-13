setMethod(
  f="summary",
  signature=c("MixClusresults"),
  function(object, ...){
      cat("**************************************************************\n\n")
      cat("DATA SET description:","\n")
      cat("* Number of individuals: ",object@data$nrow,"\n")
      cat("* Number of variables (continuous, integer, binary, ordinal): ",paste(sum(object@data$kind==1),sum(object@data$kind==2),sum(object@data$kind==3),sum(object@data$kind==4),sep=", "),"\n\n")
      cat("**************************************************************\n\n")
      cat("MODEL description:","\n")
      cat("* Number of classes:   ", length(object@param@proportions),"\n")
      cat("* Model:            ", object@model,"\n")
      cat("* Log-likelihood:  ", object@criteria$loglike,"\n")
      cat("* BIC:             ", object@criteria$bic,"\n")
      cat("* ICL:             ", object@criteria$icl,"\n")
      
      cat("**************************************************************\n\n")
      cat("PARAMETERS description:","\n")
      for (k in 1:length(object@param@proportions)){
        cat("* Class", k,"\n")         
        cat("* Proportion            ", object@param@proportions[k],"\n")
        if (any(object@data$kind==1)){
          cat("* Margins parameters of the continuous variables (Gaussian distribution):\n")
          tmp <- object@param@margins[[k]][which(object@data$kind==1),1:2]
          if (sum(object@data$kind==1)==1){
            tmp <- matrix(object@param@margins[[k]][which(object@data$kind==1),1:2],1,2)
            rownames(tmp)=colnames(object@data$data)[which(object@data$kind==1)]
          }
          colnames(tmp) <- c("Mean","Sd")
          print(tmp, digits=2)
          cat("\n")
        }
        if (any(object@data$kind==2)){
          cat("* Margins parameters of the integer variables (Poisson distribution):\n")
          tmp <- as.matrix(object@param@margins[[k]][which(object@data$kind==2),1])
          colnames(tmp) <- c("Mean")
          rownames(tmp)=rownames(object@param@correlations[[k]])[which(object@data$kind==2)]
          print(tmp, digits=2)
          cat("\n")
        }
        if (any(object@data$kind==3)){
          cat("* Margins parameters of the binary variables (Beta distribution):\n")
          tmp <- object@param@margins[[k]][which(object@data$kind==3),1]
          tmp <- cbind(tmp, 1-tmp)
          rownames(tmp) <- rownames(object@param@correlations[[k]])[which(object@data$kind==3)]
          colnames(tmp) <- c("Prob.0","Prob.1")
          print(tmp, digits=2)
          cat("\n")
        }
        if (any(object@data$kind==4)){
          cat("* Margins parameters of the ordinal variables (multinomial distribution):\n")
          tmp <- object@param@margins[[k]][which(object@data$kind==1),1:ncol(object@param@margins[[k]])]
          print(tmp, digits=2)
          cat("\n")
        }
        cat("\n")
        if (object@model=="free"){
          cat("* Correlations matrix:\n")
          print(object@param@correlations[[k]],digits=2)
        }
        cat("\n")
        
        cat("*****************************\n")
          
      }
      if (object@model!="free"){
        cat("* Correlations matrix:\n")
        print(object@param@correlations[[k]],digits=2)
        
      }

    return(invisible())
  }
)