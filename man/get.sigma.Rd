\name{get.sigma}
\alias{get.sigma}
\title{Standard deviation estimation.}
\description{This function estimates the standard deviation when \code{family="gaussian"}. This function is for internal use called by the \code{icmm} function.}
\usage{
get.sigma(Y, X, beta, alpha)
}
\arguments{
  \item{Y}{a (n*1) numeric matrix of responses.}
  \item{X}{a (n*p) numeric design matrix.}
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{alpha}{a scalar value of hyperparmeter \code{alpha}.}
}
\details{Estimate standard deviation as the mode of its full conditional distribution function when specify \code{family="gaussian"}. This function is for internal use called by the \code{icmm} function.}
\value{Return a scalar value of standard deviation.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
alpha<-0.5
# Obtain initial values from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
# Obtain sigma
sigma<-get.sigma(Y=Y, X=X, beta=beta, alpha=alpha)
}
