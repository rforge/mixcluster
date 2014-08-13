\name{MixClusanalysis}
\alias{MixClusanalysis}
\title{MixClusanalysis function}
\usage{
  MixClusanalysis(x, g, kind=NULL, model="free", nbchains=1,
                  burn=100, iterGibbs=1000, param=NULL, latent=NULL)
}
\arguments{
  \item{x}{Input data as matrix or data-frame. The binary and the ordinal variables have to use a numeric coding as follows: 0,1,2,...,number of modalities.}

  \item{g}{Integer specifying the number of classes.}
  
  \item{kind}{Vector indicating the nature of the variables as follows: 1=continuous, 2=integer, 3=binary, 4=ordinal. If this input is not specified, the function automatically detects the nature of each variables.}

  \item{model}{A vector of the competing models: "free"= mixture of Gaussian copulas without contraint, "homo"= mixture of Gaussian copulas with equal correlation matrices, "indpt"=locally independent mixture model, "all"= the three previous models.}
  
  \item{nbchains}{Number of MCMC chains.}
  
  \item{burn}{Number of iteration for the burn-in of the Gibbs sampler.}
  
  \item{iterGibbs}{Number of iteration for the parameter estimation performed via the Gibbs sampler.}
  
  \item{param}{An instance of \code{\linkS4class{MixClusparam}} to initialize the Gibbs sampler. If it is not specified, the Gibbs sampler is initialized in the maximum likelihood estimates of the locally independent mixture model.}
  
  \item{latent}{A matrix of the latent vector which can be used to initialize the Gibbs sampler. The first column indicates the class membership by using a numeric coding. The other columns indicate the values of the continuous latent variable of the Gaussian copula conditionally on the class membership.}
}
\value{
  Return an instance of \code{\linkS4class{MixClusresults}} class or a list of objects of  \code{\linkS4class{MixClusresults}} class.
}
\description{
  This function performs the cluster analysis of mixed data (continuous, integer, binary and ordinal) by using the mixture model of Gaussian copulas.
}
\examples{
# Toy example
data("simu")

mixres <- MixClusanalysis(simu, 2, model="homo", burn=50, iterGibbs=200)

## Cluster analysis of the wines data set by the three bi-components mixture models
#  data("wines")
#  mixres <- MixClusanalysis(wines, 2, model="all")

## Cluster analysis of the forest fires data set by the three components mixture model of Gaussian copulas
## without constraint on the correlation matrices
#  data(forestfires)
#  mixres <- MixClusanalysis(forestfires, 2)


}

