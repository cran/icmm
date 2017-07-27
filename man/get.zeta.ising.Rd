\name{get.zeta.ising}
\alias{get.zeta.ising}
\title{Local posterior probability estimation.}
\description{This function estimates the local posterior probability when assuming Ising prior on structured predictors.}
\usage{
get.zeta.ising(SS, beta, alpha, hyperparam, structure, edgeind, j)
}
\arguments{
  \item{SS}{a scalar value of sufficient statistic for regression coefficient.}
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{alpha}{a scalar value of hyperparameter \code{alpha}.}
  \item{hyperparam}{a two-dimensional vector of hyperparameters \code{a} and \code{b}.}
  \item{structure}{a data frame stores the information of structure among predictors.}
  \item{edgeind}{a vector stores primary keys of \code{structure}.}
  \item{j}{an index ranges from 1 to p. This function estimate a local posterior probability of predictor \code{j}.}
}
\details{Given all other parameters, this function estimates the local posterior probability or the probability that a regression coefficient is not zero conditional on other par
ameters. This function is called when assuming Ising prior on structured predictors.}
\value{Return a scalar value of local posterior probability.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
data(linearrelation)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
n<-dim(X)[1]
# Obtain initial values from lasso
data(initbetaGaussian)
initbeta<-as.matrix(initbetaGaussian)
# Get final output from ebvs
output<-icmm(Y, X, b0.start=0, b.start=initbeta, family = "gaussian", 
        ising.prior = TRUE, structure=linearrelation, estalpha = FALSE, 
        alpha = 0.5, maxiter = 100)
b0<-output$coef[1]
beta<-matrix(output$coef[-1], ncol=1)
# Get all parameters for function arguments
w<-get.wprior(beta)
alpha<-0.5
sigma<-get.sigma(Y,X,beta,alpha)
edgeind<-sort(unique(linearrelation[,1]))
hyperparam<-get.ab(beta, linearrelation, edgeind)
# Estimate local posterior probability
j<-1
Yres<-Y-b0-X\%*\%beta+X[,j]*beta[j,1]
sxy<-t(Yres)\%*\%X[,j]
ssx<-sum(X[,j]^2)
SS<-sqrt(n-1)*sxy/(sigma*ssx)
zeta<-get.zeta.ising(SS=SS, beta=beta, alpha=alpha, hyperparam=hyperparam, 
      structure=linearrelation, edgeind=edgeind, j=j)
}
