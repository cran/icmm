get.wpost<-function(SS, beta, alpha, hyperparam, structure, edgeind, j)
{
  delta<-ifelse(beta!=0,1,0)
  if (sum(edgeind==j)==0)
  {
    prob1<-(0.5*alpha*(1-pnorm(SS+alpha))*exp(alpha*SS+0.5*alpha^2)+0.5*alpha*pnorm(SS-alpha)*exp(-alpha*SS+0.5*alpha^2))*exp(-hyperparam[1])
  }
  else
  {
    Nind<-as.numeric(unlist(strsplit(structure[which(structure[,1]==j),2], ";")))
    prob1<-(0.5*alpha*(1-pnorm(SS+alpha))*exp(alpha*SS+0.5*alpha^2)+0.5*alpha*pnorm(SS-alpha)*exp(-alpha*SS+0.5*alpha^2))*exp(-hyperparam[1]-hyperparam[2]*sum(delta[Nind,1]))
  }
  prob0<-dnorm(SS)
  if (is.nan(prob1)==TRUE)
  {
    prob1<-0
  }
  wpost<-prob1/(prob0+prob1)
  return(wpost)
}
