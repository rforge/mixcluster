\name{MixClusResults-class}
\alias{MixClusResults-class}
\title{MixClusResults objects}
\description{
  This class contains all the output of the \link{MixClusClustering} function. 
  \describe{ 
    \item{model}{A list where its first element is an instance of \code{\linkS4class{MixClusModel_hetero}}, \code{\linkS4class{MixClusModel_homo}} or \code{\linkS4class{MixClusModel_indpt}}.}
    \item{param}{An instance of \code{\linkS4class{MixClusParam}} presenting the estimates.}
    \item{data}{An instance of \code{\linkS4class{MixClusData}}.}
    \item{priors}{A list of the hyper-parameters.}
    \item{initparam}{An instance of \code{\linkS4class{MixClusParam}} presenting the initial parameters.}
  }
}
\examples{
new("MixClusResults")

getSlots("MixClusResults")

}