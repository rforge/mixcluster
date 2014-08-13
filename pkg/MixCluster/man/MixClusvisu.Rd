\name{MixClusvisu}
\alias{MixClusvisu}
\title{MixClusvisu function}
\usage{
  MixClusvisu(obj, class, axe=c(1,2))
}
\arguments{
  \item{obj}{An instance of \code{\linkS4class{MixClusresults}} class}

  \item{class}{Integer specifying the class from which the PCA is performed.}
  
  \item{axe}{Factorial axes drawn by the function}
}
\value{
  This function returns an updated instance of \code{\linkS4class{MixClusresults}} class.
}
\description{
  This function updates the instance of \code{\linkS4class{MixClusresults}} class by computing the coordinates of the individuals
  on a factorial space per class (PCA). Moreover, it provides a scatterplot of the individuals and the correlation circle.
}
\examples{
data("simu")
mixres <- MixClusanalysis(simu, 2, model="homo", burn=50, iterGibbs=200)

print(mixres@criteria)
summary(mixres)

# Plot when the factorial coordinates have to be computed
mixres <- MixClusvisu(mixres, 1, c(1,2))

# Plot when the factorial coordinates are computed
mixres <- MixClusvisu(mixres, 1, c(1,3))
}

