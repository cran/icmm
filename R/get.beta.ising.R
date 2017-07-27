get.beta.ising<-function(SS, wpost, alpha, scaledfactor)
{
  positive<-1
  if (SS<0)
  {
    SS<- -SS
    positive<-0
  }

  F0<-pnorm(SS-alpha)/((1-pnorm(SS+alpha))*exp(2*alpha*SS)+pnorm(SS-alpha))
  if (wpost*F0 <= 0.5)
  {
    beta.scaled<-0
  }
  else
  {
    med1<-1-((1-pnorm(SS+alpha))*exp(2*SS*alpha)+pnorm(SS-alpha))/(2*wpost)
    med2<-SS-alpha+qnorm(med1)
    beta.scaled<-max(0, med2)
  }

  if (positive==0)
  {
    beta.scaled<- -beta.scaled
  }
  beta<-beta.scaled*scaledfactor
  return(beta)
}
