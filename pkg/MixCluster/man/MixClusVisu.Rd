\name{MixClusVisu}
\alias{MixClusVisu}
\title{MixClusVisu function}
\usage{
  MixClusVisu(obj, class, axe=c(1,2), figure, ...)
}
\arguments{
  \item{obj}{An instance of \code{\linkS4class{MixClusResults}} class}

  \item{class}{Integer specifying the class from which the PCA is performed.}
  
  \item{axe}{Factorial axes drawn by the function}
  
  \item{figure}{Which figure is plotted: "scatter"=scatterplot of the individuals, "circle"=correlation circle}
  
  \item{...}{Other inputs parameters of a function plot()}
}
\value{
 NULL.
}
\description{
  This function provides a scatterplot of the individuals on the axes obtained by a PCA performed of the 
  data set "expectation of y for a fix clas"   and the correlation circle.
}
\examples{
data("simu")

mixres <- MixClusClustering(simu, 2, nbalgo=1, burn_in=30, nbiter=100)

MixClusVisu(mixres, 1, axe=c(1,2))
}

