
\name{MixClusPrior_integer-class}
\alias{MixClusPrior_integer-class}
\title{MixClusPrior_integer objects}
\description{
  This class contains the hyper-parameters related to a Poisson distribution
  \describe{ 
    \item{ao}{First parameter of the gamma distribution.}
    \item{Ao}{Second parameter of the gamma distribution.}
  }
}
\examples{
  new("MixClusPrior_integer")
  
  getSlots("MixClusPrior_integer")
  
}