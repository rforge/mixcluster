
\name{MixClusPrior_continuous-class}
\alias{MixClusPrior_continuous-class}
\title{MixClusPrior_continuous objects}
\description{
  This class contains the hyper-parameters related to a Gaussian distribution
  \describe{ 
    \item{c0}{first parameter of the gamma distribution (related to sd)}
    \item{Co}{second parameter of the gamma distribution (related to sd)} 
    \item{bo}{first parameter of the Gaussian distribution (related to mu)}
    \item{No}{second parameter of the Gaussian distribution (related to mu)} 
  }
}
\examples{
new("MixClusPrior_continuous")

getSlots("MixClusPrior_continuous")

}