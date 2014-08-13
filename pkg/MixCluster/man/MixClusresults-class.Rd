\name{MixClusresults-class}
\alias{MixClusresults-class}
\title{MixClusresults objects}
\description{
  This class contains all the output of the MixClusanalysis function 
  \describe{ 
    \item{param}{An instance of \code{\linkS4class{MixClusparam}} class containing the parameters.}
    \item{tik}{A matrix given the posterior probabilities of the class memberships.}
    \item{partition}{A vector given the hard partition computed with the Maximum A Posteriori principle.}
    \item{criteria}{A list given the values of the following indices: the log-likelihood, BIC, ICL.}
    \item{data}{A list describing the data set.}
    \item{current_latent}{A matrix given the current values of the latent parameters.}
    \item{model}{Specification of the model (free, homoscedastic, locally independent)}
    \item{expect_y}{A list of the conditional expectation of y per class}
    \item{test_expect_y}{A boolean which indicates if the expectation of y has been computed for each class}
     
  }
}
\examples{
new("MixClusresults")

getSlots("MixClusresults")

}