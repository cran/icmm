\name{get.wprior}
\alias{get.wprior}
\title{Mixing weight estimation.}
\description{Given other parameters, this function estimates a mixing weight from the mode of its full conditional distribution function.}
\usage{
get.wprior(beta)
}
\arguments{
  \item{beta}{a (p*1) matrix of regression coefficients.}
}
\details{Given other parameters, this function estimates a mixing weight from the mode of its full conditional distribution function. This function is called when use the independent prior of predictors (no prior on structured predictors).}
\value{Return a scalar value of a mixing weight.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
# Obtain initial values from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
# Estimate the mixing weight
w<-get.wprior(beta)
}
