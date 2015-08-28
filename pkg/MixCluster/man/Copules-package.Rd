\name{MixCluster-package}
\alias{MixCluster-package}
\alias{MixCluster}
\docType{package}
\title{
MIXed-type of data CLUSTER analysis by a mixture model of Gaussian copulas
}
\description{
MixCluster performs the cluster analysis, the semi-supervised and the supervized analysis of mixed data with missing values by using a mixture of Gaussian copulas (Marbac and al, 2014).
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

The main two functions of MixCluster are \link{MixClusClustering} and \link{MixClusClassification}. They perform the parameter estimation by a Gibbs sampler then they respectively perform the cluster analysis and the semi-supervised and supervized classification.  The one-dimensional margins per class are plotted by the function \link{plot}. A summary of the results is provided by the function \link{summary}. Finally, a scatterplot of the individual using a PCA per class can be drawn by the function \link{MixClusVisu}. This function also summarizes the intra-class dependencies by drawing the correlation circle per class. Finally, MixCluster contains three data sets: one simulated data set and two real data sets used in the article Marbac and al (2014).

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

 simu[runif(300)<0.2] <- NA
# Cluster analysis by the bi-component mixture model of Gaussian copulas
# without constrain between the correlation matrices
resclustering <- MixClusClustering(simu, 2, nbalgo=1, burn_in=30, nbiter=100)

# Univariate probability distribution functions of the one-dimensional 
# margins per component weigthed by the component proportions
plot(resclustering )

# PCA-type visualization per class by using the model parameters 
# (here visualization throughout the parameters of class 1)
MixClusVisu(resclustering , 1)

# Supervized classification performed by the bi-component homoscedastic model with missing values
who <- which(runif(100)<0.6)
z <- rep(2, length(who))
z[which(who<=65)] <- 1
ressupervized  <- MixClusClassification(simu[who,], z, nbalgo=1, burn_in=30, nbiter=100)
predictsupervized <- MixClusPredict(simu[-who,], ressupervized)


# Semi-supervized classification performed by the bi-component homoscedastic model with missing values
z <- c(rep(1,65), rep(2,35))
z[-who] <- NA
ressemisupervized  <- MixClusClassification(simu, z, nbalgo=1, burn_in=30, nbiter=100)

                   
# Summary of the results
summary(ressemisupervized)



}
