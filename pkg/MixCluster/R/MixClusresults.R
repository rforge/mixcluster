setClass(
  Class="MixClusresults",
  representation=representation(
    param="MixClusparam",
    tik="matrix",
    partition="numeric",
    criteria="list",
    data="list",
    current_latent="matrix",
    model="character",
    expect_y="list",
    test_expect_y="numeric"
  ),
  prototype = prototype(
    param=GCMMparam(0,list(),list()),
    tik=matrix(nrow=0,ncol=0),
    partition=numeric(0),
    criteria=list(),
    data=list(),
    current_latent=matrix(nrow=0,ncol=0),
    model=character(0),
    expect_y=list(),
    test_expect_y=numeric(0)
  )
)


MixClusresultsctr <- function(param, tik, partition, criteria, data,  current_latent, model, expect_y, test_expect_y){
  new("MixClusresults", param=param, tik=tik, partition=partition, criteria=criteria, data=data, 
      current_latent=current_latent, model=model, expect_y=expect_y, test_expect_y=test_expect_y)
}
