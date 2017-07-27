get.pseudodata.cox<-function(Y, X, event, beta, time, ntime, sumevent)
{
  n<-dim(X)[1]
  epsilon<-10^(-5)
  hazardinc<-rep(0,ntime)
  cumuhazard<-rep(0,n)
  eta<-X%*%beta
  for(j in 1:ntime)
  {
    hazardinc[j]=sumevent[j]/sum(exp(eta[Y[,1]>=time[j]]))
  }
  # cumulative baseline hazard function
  for(i in 1:n)
  {
    cumuhazard[i]<-sum(hazardinc[time<=Y[i,1]])
  }

  mu<-(diag(cumuhazard))%*%(exp(eta))
  mu[mu[,1]==0,1]<-epsilon
  sigma2<-matrix(1/mu, ncol=1)
  z<-eta+(diag(sigma2[,1])%*%((matrix(event, ncol=1)-matrix(mu, ncol=1))))
  pseudodata<-list(z,sigma2)
  names(pseudodata)<-c("z", "sigma2")
  return(pseudodata)
}
