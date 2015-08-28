\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,MixClusResults-method}
\title{Summary of a MixClusresults class}
\arguments{
  \item{x}{A \code{\linkS4class{MixClusResults}} object}

  \item{...}{further arguments passed to or from other
  methods}
}
\value{
  NULL. Summaries to standard out.
}
\description{
  Provides result summaries of a MixClusresults class
}
\examples{
data("simu")

mixres <- MixClusClustering(simu, 2, nbalgo=1, burn_in=30, nbiter=100)

summary(mixres)

}
\seealso{
  \code{\link{summary}}
}

