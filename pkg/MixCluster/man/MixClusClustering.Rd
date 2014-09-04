\name{MixClusClustering}
\alias{MixClusClustering}
\title{MixClusClustering function}
\usage{
  MixClusClustering(x, g, model="hetero", kind=NULL, nbalgo=3,
                    burn_in=300, nbiter=300,  param)
}
\arguments{
  \item{x}{Input data as matrix or data-frame. The binary and the ordinal variables have to use a numeric coding as follows: 0,1,2,...,number of modalities.}

  \item{g}{Integer specifying the number of classes.}
  
  
  \item{model}{One of the following models: "hetero"= mixture of Gaussian copulas without constraint, "homo"= mixture of Gaussian copulas with equal correlation matrices, "indpt"=locally independent mixture model.}
 
  \item{kind}{Vector indicating the nature of the variables as follows: 1=continuous, 2=integer, 3=ordinal. If this input is not specified, the function automatically detects the nature of each variables.}
 
  \item{nbalgo}{Number of MCMC chains.}
  
  \item{burn_in}{Number of iterations for the burn-in of the Gibbs sampler.}
  
  \item{nbiter}{Number of iterations for the parameter estimation performed via the Gibbs sampler.}
  
  \item{param}{An instance of \code{\linkS4class{MixClusParam}} to initialize the Gibbs sampler. If it is not specified, the Gibbs sampler is initialized in the maximum likelihood estimates of the locally independent mixture model.}
  
 }
\value{
  Return an instance of \code{\linkS4class{MixClusResults}} class.
}
\description{
  This function performs the cluster analysis of mixed data sets with missing values by using the mixture model of Gaussian copulas.
}
\examples{
# Toy example
data("simu")

# Introduce missing values
simu[runif(300)<0.2] <- NA

# Cluster analysis by using the bi-component heteroscedastic model
mixres <- MixClusClustering(simu, 2, nbalgo=1, burn_in=30, nbiter=100)


## Cluster analysis of the wines data set by the three bi-components mixture models
#  data("wines")
#  mixres <- MixClusClustering(wines, 2)

## Cluster analysis of the forest fires data set by 
## the three components mixture model of Gaussian copulas
## without constraint on the correlation matrices
#  data(forestfires)
#  mixres <- MixClusClustering(forestfires, 2)


}

