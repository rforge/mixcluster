
\name{MixClusModel_indpt-class}
\alias{MixClusModel_indpt-class}
\title{MixClusModel_indpt objects}
\description{
  This class contains the elements of the locally independent mixture model.
  \describe{ 
    \item{g}{Class number}
    \item{nbparam}{Number of parameters required by the model.} 
    \item{loglike}{Loglikelihood.}
    \item{icl}{Loglikelihood.}
    \item{bic}{Loglikelihood.}
    \item{challenge}{Indicates the challenge: clustering, semi-supervized or supervized classification.}
  }
}
\examples{
  new("MixClusModel_indpt")
  
  getSlots("MixClusModel_indpt")
  
}