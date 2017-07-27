\name{get.ab}
\alias{get.ab}
\title{Hyperparameter estimation for \code{a} and \code{b}.}
\description{This function estimates the hyperparameters \code{a} and \code{b} for the Ising prior. This function is for internal use called by the \code{icmm} function.}
\usage{
get.ab(beta, structure, edgeind)
}
\arguments{
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{structure}{a data frame stores the information of structure among predictors.}
  \item{edgeind}{a vector stores primary keys of \code{structure}.}
}
\details{Estimate hyperparameters, \code{a} and \code{b}, using maximum pseudolikelihood estimators.}
\value{Return a two-dimensional vector where the fist element is \code{a} and the second element is \code{b}.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
data(linearrelation)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
# Suppose obtain beta from lasso
data(initbetaGaussian)
beta<-as.matrix(initbetaGaussian)
edgeind<-sort(unique(linearrelation[,1]))
hyperparameter<-get.ab(beta=beta, structure=linearrelation, edgeind=edgeind)
}
