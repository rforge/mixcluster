\name{MixClusClassification}
\alias{MixClusClassification}
\title{MixClusClassification function}
\usage{
  MixClusClassification(x, z, model="hetero", kind=NULL, nbalgo=3,
                        burn_in=300, nbiter=300,  param)
}
\arguments{
  \item{x}{Input data as matrix or data-frame. The binary and the ordinal variables have to use a numeric coding as follows: 0,1,2,...,number of modalities.}

  \item{z}{Vector of the class membership (integer between 1 and g for the labelled indivudals and NA for the unlabelled individuals. }
  
  
  \item{model}{One of the following models: "hetero"= mixture of Gaussian copulas without constraint, "homo"= mixture of Gaussian copulas with equal correlation matrices, "indpt"=locally independent mixture model.}
 
  \item{kind}{Vector indicating the nature of the variables as follows: 1=continuous, 2=integer, 3=ordinal. If this input is not specified, the function automatically detects the nature of each variables.}
 
  \item{nbalgo}{Number of MCMC chains.}
  
  \item{burn_in}{Number of iterations for the burn-in of the Gibbs sampler.}
  
  \item{nbiter}{Number of iterations for the parameter estimation performed via the Gibbs sampler.}
  
  \item{param}{An instance of \code{\linkS4class{MixClusParam}} to initialize the Gibbs sampler. If it is not specified, the Gibbs sampler is initialized in the maximum likelihood estimates of the locally independent mixture model.}
  
 }
\value{
  Return an instance of \code{\linkS4class{MixClusResults}} class.
}
\description{
  This function performs the semi-supervized and the supervized classification of mixed data sets with missing values by using the mixture model of Gaussian copulas.
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


# Semi-supervized classification performed by the bi-component homoscedastic model with missing values
z <- c(rep(1,65), rep(2,35))
z[-who] <- NA
ressemisupervized  <- MixClusClassification(simu, z, nbalgo=1, burn_in=30, nbiter=100)

                   
# Summary of the results
summary(ressemisupervized)


}

