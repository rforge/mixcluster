\name{MixClusPredict}
\alias{MixClusPredict}
\title{MixClusPredict function}
\usage{
  MixClusPredict(x, rule)
}
\arguments{
  \item{x}{Input data as matrix or data-frame to classify. The binary and the ordinal variables have to use a numeric coding as follows: 0,1,2,...,number of modalities.}

  \item{rule}{An instance of the  \code{\linkS4class{MixClusResults}} S4 class which determines the classification rule.}
  
 
 }
\value{
  Return an instance of \code{\linkS4class{MixClusResults}} class.
}
\description{
  This function performs the prediction of the class membership for a data set according to a specified model.
}
\examples{
data("simu")

simu[runif(300)<0.2] <- NA

# Supervized classification performed by the bi-component homoscedastic model with missing values
who <- which(runif(100)<0.6)
z <- rep(2, length(who))
z[which(who<=65)] <- 1
ressupervized  <- MixClusClassification(simu[who,], z, nbalgo=1, burn_in=30, nbiter=100)
predictsupervized <- MixClusPredict(simu[-who,], ressupervized)

}

