\name{simCox}
\alias{simCox}
\docType{data}
\title{Simulated data from Cox's regression model}
\description{Simulated data from Cox's regression model. A data frame with 100 observations and 402 variables. The included variables are \cr
\code{V1} A numeric vector of responses for right censored data. \cr
\code{V2} A numeric vector of status indicator: \code{0}=right censored, \code{1}=event at time \code{V1}. \cr
\code{V3}-\code{V402} 400 vectors of covariates.
} 
\usage{data(simCox)}
\format{
  A data frame of simulated data from Cox's regression model with 100 observations and 402 variables.}
\examples{
  data(simCox)
  Y<-as.matrix(simCox[,1])
  event<-as.matrix(simCox[,2])
  X<-as.matrix(simCox[,-(1:2)])
}
