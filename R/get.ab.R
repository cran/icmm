get.ab<-function(beta, structure, edgeind)
{
  delta<-ifelse(beta!=0,1,0)
  p<-dim(beta)[1]
  if (sum(delta)>0)
  {
    yy<-delta
    xx<-matrix(rep(0,p), ncol=1)
    for (j in as.numeric(edgeind))
    {
      Nind<-as.numeric(unlist(strsplit(structure[which(structure[,1]==j),2], ";")))
      xx[j,1]<-sum(delta[Nind,1])
    }
    tmpdata<-cbind(yy,xx)
    tmpdata<-as.data.frame(tmpdata)
    logitmodel<-glm(yy~xx, data=tmpdata, family=binomial)
    hyperparam<-rep(0,2)
    hyperparam<- -logitmodel$coefficients
    hyperparam[is.na(hyperparam)]<-0
  }
  else
  {
    hyperparam<-rep(0,2)
  }
  return(hyperparam)
}
