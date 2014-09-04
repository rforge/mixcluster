
\name{MixClusModel_homo-class}
\alias{MixClusModel_homo-class}
\title{MixClusModel_homo objects}
\description{
  This class contains the elements of the homoscedastic mixture model of Gaussian copulas
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
new("MixClusModel_homo")

getSlots("MixClusModel_homo")

}