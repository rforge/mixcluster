\name{MixCluster-package}
\alias{MixCluster-package}
\alias{MixCluster}
\docType{package}
\title{
MIXed-type of data CLUSTER analysis by a mixture model of Gaussian copulas
}
\description{
MixCluster performs the cluster analysis of mixed data by using a mixture of Gaussian copulas (Marbac and al, 2014).
Thus, it can analyze data sets composed by continuous, integer, binary or ordinal variables.
The one-dimensional margins of each component follow classical distributions (Gaussian, Poisson, multinomial) and the intra-class dependencies are modelized. A Gibbs sampler performs the Bayesian inference. Finally, a PCA-type approach per class allows to visualize the individuals per class and summarizes the intra-class dependencies.
}
\details{
\tabular{ll}{
Package: \tab MixCluster \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-11\cr
License: \tab GPL (>=2)\cr
}

The main function of MixCluster is \link{MixClusanalysis} and performs the parameter estimation by a Gibbs sampler. It requires to mandatory input parameters: a data set and a class number. Note that the data set is composed by objects (in row) which are described by different variables (in column) using a numeric coding. The one-dimensional margins per class are plotted by the function \link{plot}. A summary of the results is provided by the function \link{summary}. Finally, a scatterplot of the individual using a PCA per class can be drawn by the function \link{MixClusvisu}. This function also summarizes the intra-class dependencies by drawing the correlation circle per class. Finally, MixCluster contains three data sets: one simulated data set and two real data sets used in the article Marbac and al (2014).

}
\author{
Matthieu Marbac & Christophe Biernacki & Vincent Vandewalle

Maintainer: Matthieu Marbac <matthieu.marbac-lourdelle@inria.fr>
}
\references{
M. Marbac, C. Biernacki & V. Vandewalle (2014). Model-based clustering of Gaussian copulas for mixed data. Preprint.
}



\examples{
# Importation of the data set
data("simu")

# Cluster analysis by the bi-component mixture model of Gaussian copulas
# with equal correlation matrices between classes
mixres <- MixClusanalysis(simu, 2, model="homo", burn=50, iterGibbs=200)

# Summary of the results
summary(mixres)

# Univariate probability distribution functions of the one-dimensional 
# margins per component weigthed by the component proportions
plot(mixres)

# PCA-type visualization per class by using the model parameters 
# (here visualization throughout the parameters of class 1)
mixres <- MixClusvisu(mixres, 1)

}
