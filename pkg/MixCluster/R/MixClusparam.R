setClass(
  Class="MixClusparam",
  representation=representation(
    proportions="numeric",
    margins="list",
    correlations="list"
  ),
  prototype = prototype(
    proportions=numeric(0),
    margins=list(),
    correlations=list()
  )
)


GCMMparam <- function(proportions,margins,correlations){
  new("MixClusparam",proportions=proportions,margins=margins,correlations=correlations)
}
