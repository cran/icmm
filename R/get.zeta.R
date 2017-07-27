get.zeta<-function(SS, w, alpha)
{
  betaSS<-beta.laplace(ifelse(SS>=38, 38, SS), a=alpha)
  zeta<-(1+betaSS)/((1/w)+betaSS)
  return(zeta)
}
