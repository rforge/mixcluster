
\name{MixClusModel_hetero-class}
\alias{MixClusModel_hetero-class}
\title{MixClusModel_hetero objects}
\description{
  This class contains the elements of the heteroscedastic mixture model of Gaussian copulas
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
new("MixClusModel_hetero")

getSlots("MixClusModel_hetero")

}