\name{get.wpost}
\alias{get.wpost}
\title{Estimate posterior probability of mixing weight.}
\description{With the Ising prior on structured predictors, this function gets the posterior probability of mixing weight.}
\usage{
get.wpost(SS, beta, alpha, hyperparam, structure, edgeind, j)
}
\arguments{
  \item{SS}{a scalar value of sufficient statistic for regression coefficient.}
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{alpha}{a scalar value of hyperparameter \code{alpha}.}
  \item{hyperparam}{a two-dimensional vector of hyperparameters \code{a} and \code{b}.}
  \item{structure}{a data frame stores the information of structure among predictors.}
  \item{edgeind}{a vector stores primary keys of \code{structure}.}
  \item{j}{an index ranges from 1 to p. This function estimates a  posterior probability of a mixing weight corresponding to predictor \code{j}.}
}
\details{With the Ising prior on structured predictors, the problem is transformed into the realm of empirical Bayes thresholding with Laplace prior by estimating the posterior probability of mixing weight. The posterior probability is used to find the posterior median of a regression coefficient.}
\value{Return a scalar value of a posterior probability of mixing weight for predictor.}
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
\donttest{sigma<-get.sigma(Y=Y, X=X, beta=beta, alpha=alpha)}
\donttest{hyperparam<-get.ab(beta, linearrelation, edgeind)}
# Estimate the posterior probability of first predictor
j<-1
Yres<-Y-X\%*\%beta+X[,j]*beta[j,1]
sxy<-t(Yres)\%*\%X[,j]
ssx<-sum(X[,j]^2)
\donttest{SS<-sqrt(n-1)*sxy/(sigma*ssx)}
\donttest{wpost<-get.wpost(SS=SS, beta=beta, alpha=alpha, hyperparam=hyperparam, 
       structure=linearrelation, edgeind=edgeind, j=j)}
}
