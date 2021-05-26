\name{get.zeta}
\alias{get.zeta}
\title{Local posterior probability estimation}
\description{This function estimates the local posterior probability when assuming no prior on structured predictors.}
\usage{
get.zeta(SS, w, alpha)
}
\arguments{
  \item{SS}{a scalar value of sufficient statistic for regression coefficient.}
  \item{w}{a scalar value of mixing weight.}
  \item{alpha}{a scalar value of hyperparameter \code{alpha}.}
}
\details{Given all other parameters, this function estimates the local posterior probability or the probability that a regression coefficient is not zero conditional on other parameters. This function is called when assuming no prior on structured predictors.}
\value{Return a scalar value of local posterior probability.}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\examples{
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
n<-dim(X)[1]
# Obtain initial values from lasso
data(initbetaGaussian)
initbeta<-as.matrix(initbetaGaussian)
# Obtain the final output from ebvs
\donttest{output<-icmm(Y, X, b0.start=0, b.start=initbeta, family = "gaussian", 
        ising.prior = FALSE, estalpha = FALSE, alpha = 0.5, maxiter = 100)}
\donttest{b0<-output$coef[1]}
\donttest{beta<-matrix(output$coef[-1], ncol=1)}
# Get all parameters for function arguments
\donttest{w<-get.wprior(beta)}
\donttest{alpha<-0.5}
\donttest{sigma<-get.sigma(Y,X,beta,alpha)}
# Estimate local posterior probability
\donttest{j<-1}
\donttest{Yres<-Y-b0-X\%*\%beta+X[,j]*beta[j,1]}
\donttest{sxy<-t(Yres)\%*\%X[,j]}
\donttest{ssx<-sum(X[,j]^2)}
\donttest{SS<-sqrt(n-1)*sxy/(sigma*ssx)}
\donttest{zeta<-get.zeta(SS=SS, w=w, alpha=alpha)}
}
