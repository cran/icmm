\name{get.beta}
\alias{get.beta}
\title{Obtain model coefficient without assuming prior on structure of predictors.}
\description{Given a sufficient statistic for a regression coefficient, this funciton estimates a regression coefficient without assuming prior on structure of predictors.}
\usage{
get.beta(SS, w, alpha, scaledfactor)
}
\arguments{
  \item{SS}{a scalar value of sufficient statistic for a regression coefficient.}
  \item{w}{a scalar value of mixing weight.}
  \item{alpha}{a scalar value of hyperparameter \code{alpha}.}
  \item{scaledfactor}{a scalar value of multiplicative factor.}
}
\details{Empirical Bayes thresholding is employed to obtain a posterior median of a regression coefficient.}
\value{a scalar value of regression coefficient.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
n<-dim(X)[1]
# Obtain initial values from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
# Initiate all other parameters
w<-0.5
alpha<-0.5
sigma<-get.sigma(Y=Y, X=X, beta=beta, alpha=alpha)
# Obtain a sufficient statistic
j<-1
Yres<-Y-X\%*\%beta+X[,j]*beta[j,1]
sxy<-t(Yres)\%*\%X[,j]
ssx<-sum(X[,j]^2)
SS<-sqrt(n-1)*sxy/(sigma*ssx)
beta[j,1]<-get.beta(SS=SS, w=w, alpha=alpha, scaledfactor=sigma/sqrt(n-1))
}
