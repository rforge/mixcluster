\name{MixCluster-package}
\alias{MixCluster-package}
\alias{MixCluster}
\docType{package}
\title{
Mixed data Clustering by a mixture model of Gaussian copulas
}
\description{
MixCluster performs the cluster analysis of mixed data with missing values by using a mixture of Gaussian copulas (Marbac and al, 2015).
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

The main two functions of MixCluster is \link{MixClusClustering}. It performs the parameter estimation by a Gibbs sampler then they respectively perform the cluster analysis. The one-dimensional margins per class are plotted by the function \link{plot}. A summary of the results is provided by the function \link{summary}. Finally, a scatterplot of the individual using a PCA per class can be drawn by the function \link{MixClusVisu}. This function also summarizes the intra-class dependencies by drawing the correlation circle per class. Finally, MixCluster contains three data sets: one simulated data set and two real data sets used in the article Marbac and al (2015).
}
\author{
Matthieu Marbac & Christophe Biernacki & Vincent Vandewalle

Maintainer: Matthieu Marbac <matthieu.marbac@gmail.com>
}
\references{
M. Marbac, C. Biernacki & V. Vandewalle (2015). Model-based clustering of Gaussian copulas for mixed data. Preprint.
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
