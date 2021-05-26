\name{get.pseudodata.binomial}
\alias{get.pseudodata.binomial}
\title{Obtain pseudodata based on the binary logistic regression model.}
\description{For logistic regression, given the current estimates of regression coefficients, working responses and their corresponding weights are obtained.}
\usage{
get.pseudodata.binomial(Y, X, beta0, beta, niter)
}
\arguments{
  \item{Y}{an (n*1) numeric matrix of responses.}
  \item{X}{an (n*p) numeric design matrix.}
  \item{beta0}{a scalar value of intercept term.}
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{niter}{number of iterations in ICM/M algorithm.}
}
\value{Return a list including elements
  \item{z}{an (n*1) matrix of working responses}
  \item{sigma2}{an (n*1) matrix of inverse of weights.}
}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simBinomial)
Y<-as.matrix(simBinomial[,1])
X<-as.matrix(simBinomial[,-1])
p<-dim(X)[2]
# Obtain initial values from lasso
data(initbetaBinomial)
initbeta<-as.matrix(initbetaBinomial)
# Get Pseudodata
\donttest{pseudodata<-get.pseudodata.binomial(Y=Y, X=X, beta0=0, beta=initbeta, niter=1)}
\donttest{z<-pseudodata$z}
\donttest{sigma<-sqrt(pseudodata$sigma2)}
}
