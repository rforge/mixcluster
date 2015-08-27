\name{MixClusparam-class}
\alias{MixClusparam-class}
\title{MixClusparam objects}
\description{
  This class contains the parameters of the mixture model of Gaussian copulas
  function 
  \describe{ 
    \item{proportions}{A numeric vector indicating the class proportions.}
    \item{margins}{A list indicating the parameters of the one-dimensional marings per class. Parameters of one class are displayed in a matrix whose the parameters of the variable j are presented in the row j.}
    \item{correlations}{A list indicating the correlation matrices per class.} 
  }
}
\examples{
new("MixClusparam")

getSlots("MixClusparam")

}