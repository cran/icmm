get.alpha<-function(beta, scaledfactor)
{
  nz<-sum(beta!=0)
  alpha<-nz*scaledfactor
  if (is.nan(alpha)==TRUE)
  {
    alpha<-0.5
  }
  return(alpha)
}
