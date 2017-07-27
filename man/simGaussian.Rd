\name{simGaussian}
\alias{simGaussian}
\docType{data}
\title{Simulated data from the normal linear regression model}
\description{Simulated data from the normal linear regression model. A data frame with 100 observations and 401 variables. The included variables are\cr
\code{V1} A numeric vector of responses.\cr
\code{V2-V401} 400 vectors of covariates.
} 
\usage{data(simGaussian)}
\format{
  A data frame of simulated data from the normal linear regression with 100 observations and 401 variables.
}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
}

