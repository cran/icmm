get.wprior<-function(beta)
{
  p<-dim(beta)[1]
  w<-sum(beta!=0)/p
  return(w)
}
