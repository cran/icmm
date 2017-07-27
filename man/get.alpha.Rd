\name{get.alpha}
\alias{get.alpha}
\title{Hyperparameter estimation for \code{alpha}.}
\description{This function estimates a hyperparameter \code{alpha}, a scale parameter in Laplace denstiy. This function is for internal use called by the \code{icmm} function.}
\usage{
get.alpha(beta, scaledfactor)
}
\arguments{
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{scaledfactor}{a scalar value of multiplicative factor.}
}
\details{This function estimates a hyperparameter \code{alpha}, a scale parameter in Laplace density as the mode of its full conditional distribution function.}
\value{Return a scalar value of \code{alpha}.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
n<-dim(X)[1]
# Obtain initial values of beta from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
# Initiate alpha
alpha<-0.5
# Estimate sigma
e<-Y-X\%*\%beta
nz<-sum(beta[,1]!=0)
sigma<-get.sigma(Y=Y, X=X, beta=beta, alpha=alpha)
# Update alpha as the mode of its full conditional distribution function
alpha<-get.alpha(beta=beta, scaledfactor=1/(sqrt(n-1)*sum(abs(beta))/sigma))
}
