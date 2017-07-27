get.beta<-function(SS, w, alpha, scaledfactor)
{
  beta.scaled<-postmed(SS,w=w, prior="laplace", a=alpha)
  beta<-beta.scaled*scaledfactor
  return(beta)
}
