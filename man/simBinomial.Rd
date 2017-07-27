\name{simBinomial}
\alias{simBinomial}
\docType{data}
\title{Simulated data from the binary logistic regression model}
\description{Simulated data from the binary logistic regression model. A data frame with 100 observations and 401 variables. The included variables are \cr
\code{V1} A numeric vector of binary responses where each entry is either \code{0} or \code{1}. \cr
\code{V2-V401} 400 vectors of covariates.
} 
\usage{data(simBinomial)}
\format{
  A data frame of simulated data from the binary logistic regression with 100 observations and 401 variables.
}
\examples{
  data(simBinomial)
  Y<-as.matrix(simBinomial[,1])
  X<-as.matrix(simBinomial[,-1])
}
