\name{MixClusData-class}
\alias{MixClusData-class}
\title{MixClusData objects}
\description{
  This class contains the elements related to the data set.
  \describe{ 
    \item{x}{A matrix presenting the observed data. Missing values are indicated by NA.}
    \item{o}{A list where element j contains the individual whose variable j is observed.}
    \item{e}{Number of observed variables.} 
    \item{n}{Number of individuals.} 
    \item{kind}{A vector containing the nature of the variables (1: continuous, 2: integer, 3: ordinal).} 
    \item{partition}{A vector indicating the partition estimated by applying the Maximum A posteriori rule.} 
    \item{tik}{A matrix where row i and column k denotes the probability that individual i arises from class k.} 
    \item{proba}{A matrix where row i and column k denotes the pdf of individual i under component k.} 
    \item{condexpec}{A list containing the expectation of y per class.} 
  }
}
\examples{
new("MixClusData")

getSlots("MixClusData")

}
