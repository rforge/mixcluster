\name{MixClusUpdateForVisu}
\alias{MixClusUpdateForVisu}
\title{MixClusUpdateForVisu function}
\usage{
  MixClusUpdateForVisu(output)
}
\arguments{
  \item{output}{An instance of \code{\linkS4class{MixClusResults}} class}
}
\value{
 An instance of \code{\linkS4class{MixClusResults}} class
}
\description{
  This function provides computes the conditional expectation of the latent vectors related to the Gaussian copulas. It must be used before calling the function \link{MixClusVisu}.
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

