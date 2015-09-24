\name{MixClusClustering}
\alias{MixClusClustering}
\title{MixClusClustering function}
\usage{
  MixClusClustering(x, g, model="hetero", kind=NULL, nbalgo=1, burn_in=100, nbiter=1000, param)
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
\dontrun{
# Loading of a dataset simulated from a bi-component mixture model of Gaussian copulas
# (see Example 2.2 page 6)
# The first column indicates the class membership
# The last three column are used for the clustering
data(simu)

# Cluster analysis by the bi-component mixture model of Gaussian copulas
# without constrain between the correlation matrices
res.mixclus <- MixClusClustering(simu[,-1], 2)

# Confusion matrix between the estimated (row) and the true (column) partition
table(res.mixclus@data@partition, simu[,1])

# Summary of the model
summary(res.mixclus)

# Visualisation
# Update of the results (computing the conditional expectations of the latent vectors
# related to the Gaussian copulas)
res.mixclus <- MixClusUpdateForVisu(res.mixclus)

# Scatterplot of the individuals  (Figure 1.(a)) described by three variables:
# one continuous (abscissa), one integer (ordiate) and one binary (symbol).
# Colors indicate the component memberships
plot(simu[,2:3], col=simu[,1], pch=16+simu[,4], xlab=expression(x^1), ylab=expression(x^2))

# Scatterplot of the individuals in the first PCA-map of the first-component of the model
MixClusVisu(res.mixclus, class = 2, figure = "scatter", xlim=c(-10,4), ylim=c(-4,4))

# Correlation circle of the first PCA-map of the first-component of the model
MixClusVisu(res.mixclus, class = 2, figure = "circle")

}
}

