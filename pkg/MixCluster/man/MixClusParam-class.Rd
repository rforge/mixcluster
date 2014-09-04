\name{MixClusParam-class}
\alias{MixClusParam-class}
\title{MixClusParam objects}
\description{
  This class contains the parameters of the mixture model of Gaussian copulas
  function 
  \describe{ 
    \item{pi}{A numeric vector indicating the class proportions.}
    \item{beta}{A list indicating the parameters of the one-dimensional marings per class. Parameters of one class are displayed in a list whose element j groups the parameters related to variable j and are an instance of \code{\linkS4class{MixClusParam_continuous}}, \code{\linkS4class{MixClusParam_integer}} or \code{\linkS4class{MixClusParam_ordinal}}.}
    \item{correlations}{A list indicating the correlation matrices per class.} 
  }
}
\examples{
new("MixClusParam")

getSlots("MixClusParam")

}
