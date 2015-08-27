\docType{methods}
\name{summary}
\alias{summary}
\alias{summary,MixClusresults-method}
\title{Summary of a MixClusresults class}
\arguments{
  \item{x}{A \code{\linkS4class{MixClusresults}} object}

  \item{...}{further arguments passed to or from other
  methods}
}
\value{
  NULL. Summaries to standard out.
}
\description{
  Produce result summaries of a MixClusresults class
}
\examples{
data("simu")
mixres <- MixClusanalysis(simu, 2, model="homo", burn=50, iterGibbs=200)

summary(mixres)

}
\seealso{
  \code{\link{summary}}
}

