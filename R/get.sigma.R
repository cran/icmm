get.sigma<-function(Y,X,beta,alpha)
{
  e<-Y-X%*%beta
  n<-dim(X)[1]
  nz<-sum(beta[,1]!=0)
  sigma<-(alpha*sqrt(n-1)*sum(abs(beta))+sqrt(alpha^2*(n-1)*(sum(abs(beta)))^2+4*(n+nz+1)*sum(e^2)) )/(2*(n+nz+1))
  return(sigma)
}
