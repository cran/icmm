\name{get.beta.ising}
\alias{get.beta.ising}
\title{Obtain a regression coefficient when assuming Ising prior (with structured predictors).}
\description{Given a sufficient statistic for a regression coefficient, this function estimates a coefficient when assuming the Ising model to incorporate the prior of structured predictors.}
\usage{
get.beta.ising(SS, wpost, alpha, scaledfactor)
}
\arguments{
  \item{SS}{a sufficient statistic for a regression coefficient.}
  \item{wpost}{a posterior probability of mixing weight.}
  \item{alpha}{a scalar value for hyperparameter \code{alpha}.}
  \item{scaledfactor}{a scalar value for multiplicative factor.}
}
\details{Given a posterior probability of mixing weight, empirical Bayes thresholding is employed to obtain a posterior median of a regression coefficient.}
\value{a scalar value of regression coefficient.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
n<-dim(X)[1]
data(linearrelation)
edgeind<-sort(unique(linearrelation[,1]))
# Obtain initial values from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
# Initiate all other parameters
alpha<-0.5
sigma<-get.sigma(Y=Y, X=X, beta=beta, alpha=alpha)
hyperparam<-get.ab(beta, linearrelation, edgeind)
# Obtain regression coefficient
j<-1
Yres<-Y-X\%*\%beta+X[,j]*beta[j,1]
sxy<-t(Yres)\%*\%X[,j]
ssx<-sum(X[,j]^2)
SS<-sqrt(n-1)*sxy/(sigma*ssx)
wpost<-get.wpost(SS, beta, alpha, hyperparam, linearrelation, edgeind, j)
beta[j,1]<-get.beta.ising(SS=SS, wpost=wpost, alpha=alpha, 
           scaledfactor=sigma/sqrt(n-1))
}
