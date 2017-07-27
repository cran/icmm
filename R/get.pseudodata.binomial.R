get.pseudodata.binomial<-function(Y, X, beta0, beta, niter)
{
  epsilon<-10^(-5)
  n<-dim(X)[1]
  eta<-beta0+X%*%beta
  eta[eta[,1]>30,1]<-30
  eta[eta[,1]< -30,1]<- -30
  prob<-exp(eta)/(1+exp(eta))
  if (niter==0)
  {
    prob[Y[,1]==0,1]<-0.25
    prob[Y[,1]==1,1]<-0.75
  }
  indzero<-which(prob[,1]<=epsilon)
  indone<-which(prob[,1]>=(1-epsilon))
  indpositive<-which(eta>0)
  indnegative<-which(eta<0)

  hazard_positive<-rep(1,n)
  hazard_positive[indpositive]<-1/(1+exp(-eta[indpositive,1]))
  hazard_positive[indnegative]<-exp(eta[indnegative,1])/(exp(eta[indnegative,1])+1)

  hazard_negative<-rep(1,n)
  hazard_negative[indpositive]<-exp(-eta[indpositive,1])/(1+exp(-eta[indpositive,1]))
  hazard_negative[indnegative]<-1/(1+exp(eta[indnegative,1]))

  prob[indzero,1]<-0
  prob[indone,1]<-1

  z<-eta+(Y-prob)/(prob*(1-prob))
  z[indone,1]<-eta[indone,1]+1/hazard_positive[indone]
  z[indzero,1]<-eta[indzero,1]-1/hazard_negative[indzero]

  sigma2<-1/(hazard_positive*hazard_negative)
  sigma2<-as.matrix(sigma2)
  z<-matrix(z, ncol=1)
  pseudodata<-list(z,sigma2)
  names(pseudodata)<-c("z", "sigma2")
  return(pseudodata)
}
