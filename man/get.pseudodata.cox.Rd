\name{get.pseudodata.cox}
\alias{get.pseudodata.cox}
\title{Obtain pseudodata based on the Cox's regression model.}
\description{For Cox's regression model, given the current estimates of regression coefficients, working responses and their corresponding weights are obtained.}
\usage{
get.pseudodata.cox(Y, X, event, beta, time, ntime, sumevent)
}
\arguments{
  \item{Y}{a (n*1) numeric matrix of time response.}
  \item{X}{a (n*p) numeric design matrix.}
  \item{event}{a (n*1) numeric matrix of status: of status indicator: \code{0}=right censored, \code{1}=event at time.}
  \item{beta}{a (p*1) matrix of regression coefficients.}
  \item{time}{a vector or sorted value of \code{Y}.}
  \item{ntime}{length of the vector \code{time}.}
  \item{sumevent}{a vector of size \code{ntime} where each element is the sum of event where \code{Y} is equal to each value in \code{time}.}
}
\value{Return a list including elements
  \item{z}{a (n*1) matrix of working responses}
  \item{sigma2}{a (n*1) matrix of inverse of weights.}
}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simCox)
Y<-as.matrix(simCox[,1])
event<-as.matrix(simCox[,2])
X<-as.matrix(simCox[,-(1:2)])
time<-sort(unique(Y))
ntime<-length(time)
# sum of event_i where y_i =time_k
sumevent<-rep(0, ntime)
for(j in 1:ntime)
{
  sumevent[j]<-sum(event[Y[,1]==time[j]])
}
# Obtain initial values from lasso
data(initbetaCox)
initbeta<-as.matrix(initbetaCox)
# Get Pseudodata
pseudodata<-get.pseudodata.cox(Y, X, event, initbeta, time, ntime, sumevent)
z<-pseudodata$z
sigma<-sqrt(pseudodata$sigma2)
}
