\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,MixClusresults-method}
\title{Plot the probability distribution functions of the one-dimensional margins per components}
\arguments{
  \item{x}{A \code{\linkS4class{MixClusresults}} object}

  \item{y}{A numeric vector indicating the plotted variables}
}
\value{
  NULL.
}
\description{
  This function plots the probability distribution functions of the one-dimensional margins per components
  weigthed by the component proportions.
}
\examples{
# Importation of the data set
data("simu")

# Cluster analysis by the bi-component mixture model of Gaussian copulas
# with equal correlation matrices between classes
mixres <- MixClusanalysis(simu, 2, model="homo", burn=50, iterGibbs=200)

# Univariate probability distribution function of the one-dimensional 
# margins per component weigthed by the component proportions
plot(mixres)


}
\seealso{
  \code{\link{plot}}
}

