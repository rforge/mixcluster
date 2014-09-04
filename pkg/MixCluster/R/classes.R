# Definition of the class related to the data
# x: data set, the ordinal variables are coding as follows {0,1...,m-1}
# o: boolean matrix 1: observed value, 0: missing value
# e: number of variables
# n: number of individuals
# kind: kind of variables
setClass(Class="MixClusData",
  representation=representation(x="matrix", o="list", e="numeric", n="numeric",
                                kind="numeric", partition="numeric", tik="matrix", proba="matrix", condexpec="list"),
  prototype = prototype(x=matrix(0,0,0), o=list(), e=numeric(0), n=numeric(0), 
                        kind=numeric(0), partition=numeric(0), tik=matrix(0,0,0), proba=matrix(0,0,0), condexpec=list() )
)




# Definition of the class related to the latent vectors (z and y)
setClass(Class="MixClusLatent",
  representation=representation(y="matrix", z="numeric"),
  prototype = prototype(y=matrix(0,0,0), z=numeric(0))
)


# Definition of the class related to the parameters
# pi: proportions
# beta: all the margin parameters
# beta[[k]]: all the margin parameters of component k
# beta[[k]][[j]]: margin parameters of margin j for component k (beta[[k]][[j]] is an instance of the following three
# classes: MixClusParam_continuous, MixClusParam_integer, MixClusParam_ordinal)
# correl: list of the correlation matrices
setClass(Class="MixClusParam",
         representation=representation(pi="numeric", beta="list", correl="list"),
         prototype = prototype(pi=numeric(0), beta=list(), correl=list())
)

setClass(Class="MixClusParam_continuous",
         representation=representation(mu="numeric", sd="numeric"),
         prototype = prototype(mu=numeric(0), sd=numeric(0))
)

setClass(Class="MixClusParam_integer",
         representation=representation(beta="numeric"),
         prototype = prototype(beta=numeric(0))
)

setClass(Class="MixClusParam_ordinal",
         representation=representation(beta="numeric"),
         prototype = prototype(beta=numeric(0))
)





MixClusPrior <- function(data, param){
  prior <- list()
  for (j in 1:data@e)
    prior[[j]] <- Findprior(data@x[data@o[[j]],j], param@beta[[1]][[j]])
  return(prior)
}



setClass(Class="MixClusPrior_continuous",
         representation=representation(co="numeric", Co="numeric", bo="numeric", No="numeric"),
         prototype = prototype( co=numeric(0), Co=numeric(0), bo=numeric(0), No=numeric(0))
)



setClass(Class="MixClusPrior_integer",
         representation=representation(ao="numeric", Ao="numeric"),
         prototype = prototype( ao=numeric(0), Ao=numeric(0))
)



setClass(Class="MixClusPrior_ordinal",
         representation=representation(a="numeric"),
         prototype = prototype( a=numeric(0))
)


setClass(Class="MixClusModel_hetero",
         representation=representation(g="numeric", nbparam="numeric", loglike="numeric", bic="numeric", icl="numeric", challenge="character"),
         prototype = prototype( g=numeric(0), nbparam=numeric(0), loglike=numeric(0), bic=numeric(0), icl=numeric(0), challenge=character(0))
)
setClass(Class="MixClusModel_homo",
         representation=representation(g="numeric", nbparam="numeric", loglike="numeric", bic="numeric", icl="numeric", challenge="character"),
         prototype = prototype( g=numeric(0), nbparam=numeric(0), loglike=numeric(0), bic=numeric(0), icl=numeric(0), challenge=character(0))
)
setClass(Class="MixClusModel_indpt",
         representation=representation(g="numeric", nbparam="numeric", loglike="numeric", bic="numeric", icl="numeric", challenge="character"),
         prototype = prototype( g=numeric(0), nbparam=numeric(0), loglike=numeric(0), bic=numeric(0), icl=numeric(0), challenge=character(0))
)



setClass(Class="MixClusResults",
         representation=representation(model="list", param="MixClusParam", data="MixClusData", priors="list", initparam="MixClusParam"),
         prototype = prototype(model=list(), param=new("MixClusParam"), data=new("MixClusData"), priors=list(), initparam=new("MixClusParam"))
)


MixClusModel <- function(model, g, data, challenge){
  if (model=="hetero"){
    nbparam <- (g-1) + g*(data@e*(data@e-1)/2 ) + 2*g*sum(data@kind==1) + g*sum(data@kind==2)
    if (any(data@kind>2)){
      for (j in which(data@kind>2))
        nbparam <- nbparam + g*(length(unique(data@x[data@o[[j]],j]))-1)
    }
    model <- new("MixClusModel_hetero", g=g, nbparam=nbparam, loglike=-Inf, bic=-Inf, icl=-Inf, challenge=challenge)
  }else if (model=="homo"){
    nbparam <- (g-1) + (data@e*(data@e-1)/2 ) + 2*g*sum(data@kind==1) + g*sum(data@kind==2)
    if (any(data@kind>2)){
      for (j in which(data@kind>2))
        nbparam <- nbparam + g*(length(unique(data@x[data@o[[j]],j]))-1)
    }
    
    model <- new("MixClusModel_homo", g=g, nbparam=nbparam, loglike=-Inf, bic=-Inf, icl=-Inf, challenge=challenge)
  }else{
    nbparam <- (g-1) + 2*g*sum(data@kind==1) + g*sum(data@kind==2)
    if (any(data@kind>2)){
      for (j in which(data@kind>2))
        nbparam <- nbparam + g*(length(unique(data@x[data@o[[j]],j]))-1)
    }
    model <- new("MixClusModel_indpt", g=g, nbparam=nbparam, loglike=-Inf, bic=-Inf, icl=-Inf, challenge=challenge)
  }
  return(model)
}

