\name{linearrelation}
\alias{linearrelation}
\docType{data}
\title{Linear structure of predictors}
\description{This data frame is used as an example to store the structure of predictors or the edge set of an undirected graph. For this data frame, the linear chain is assumed for each predictor.}
\usage{data(linearrelation)}
\format{
  A data frame with 400 observations and 2 variables as follows.
  \describe{
    \item{\code{Index}}{an index of the predictor/node which has at least one edge.}
    \item{\code{EdgeIndices}}{a string of all indices having an edge connected to \code{Index} separated by semicolon(;).}
  }
}
\details{This structure of predictors assumes a linear chain for each predictor which its immediate neighbors. For example, j-predictor is connected to (j-1)-predictor and (j+1)-predictor. The example for the entry in the data frame is \code{Index="5"} and \code{EdgeIndices="4;6"}.}
\examples{
data(linearrelation)
# To see the format of linearrelation data frame
head(linearrelation)
}
