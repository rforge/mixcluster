sample_proportions <- function(v){
  return(as.numeric(rdirichlet(1,v-0.5)))
}