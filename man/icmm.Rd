\name{icmm}
\alias{icmm}
\title{Empirical Bayes Variable Selection}
\description{Empirical Bayes variable selection via the ICM/M algorithm.}
\usage{
icmm(Y, X, event, b0.start, b.start, family = "gaussian", 
     ising.prior = FALSE, structure, estalpha = FALSE, 
     alpha = 0.5, maxiter = 100)
}
\arguments{
  \item{Y}{an (n*1) numeric matrix of responses.}
  \item{X}{an (n*p) numeric design matrix.}
  \item{event}{an (n*1) numeric matrix of status for censored data: \code{0}=censored data, \code{1}=event at time. \code{event} is required when \code{family="cox"}.}
  \item{b0.start}{a starting value of intercept term (optional).}
  \item{b.start}{a (p*1) matrix of starting values for regression coefficients.}
  \item{family}{specification of the model. It can be one of these three models: \code{"gaussian"}, \code{"binomial"}, \code{"cox"}. The default is \code{family="gaussian"}.}
  \item{ising.prior}{a logical flag for Ising prior utilization. \code{ising.prior=FALSE} assumes no prior on structure among predictors. \code{ising.prior=TRUE} indicates incorporation of Ising prior to capture structure among predictors in modeling process.}
  \item{structure}{a data frame stores the information of structured predictors (need to specify when \code{ising.prior=TRUE}).}
  \item{estalpha}{a logical flag specifying whether to obtain \code{alpha} via ICM/M algorithm.}
  \item{alpha}{a scalar value of scale parameter in Laplace density (non-zero part of prior). The default value is \code{0.5}.}
  \item{maxiter}{a maximum values of iterations for ICM/M algorithm.}
}
\details{The main function for empirical Bayes variable selection. Iterative conditional modes/medians (ICM/M) is implemented in this function. The basic problem is to estimate regression coefficients in high-dimensional data (i.e., large p small n) and we assume that most coefficients are zero. This function also allows the prior of structure of covariates to be incorporated in the model.}
\value{Return a list including elements
  \item{coef}{a vector of model coefficients. The first element is an intercept term when specifying \code{family="gaussian"} or \code{family="binomial"}.}
  \item{iterations}{number of iterations of ICM/M.}
  \item{alpha}{a scalar value of \code{alpha}.}
  \item{postprob}{a p-vector of local posterior probabilities or zeta.}
}
\author{Vitara Pungpapong, Min Zhang, Dabao Zhang}
\references{
Pungpapong, V., Zhang, M. and Zhang, D. (2015). Selecting massive variables using an iterated conditional modes/medians algorithm. Electronic Journal of Statistics. 9:1243-1266. <doi:10.1214/15-EJS1034>.\cr
Pungpapong, V., Zhang, M. and Zhang, D. (2020). Integrating Biological Knowledge Into Case-Control Analysis Through Iterated Conditional Modes/Medians Algorithm. Journal of Computational Biology. 27(7): 1171-1179. <doi:10.1089/cmb.2019.0319>.
}
\seealso{
\code{get.ab}, \code{get.alpha}, \code{get.beta}, \code{get.beta.ising}, \code{get.pseudodata.binomial}, \cr
\code{get.pseudodata.cox}, \code{get.sigma}, \code{get.wprior}, \code{get.zeta}, \code{get.zeta.ising}
}
\examples{
# Normal linear regression model
# With no prior on structure among predictors
data(simGaussian)
Y<-as.matrix(simGaussian[,1])
X<-as.matrix(simGaussian[,-1])
# Obtain initial values from lasso
data(initbetaGaussian)
initbeta<-as.matrix(initbetaGaussian)
\donttest{result<-icmm(Y=Y, X=X, b.start=initbeta, family="gaussian", 
        ising.prior=FALSE, estalpha=FALSE, alpha=0.5, maxiter=100)}
\donttest{result$coef}
\donttest{result$iterations}
\donttest{result$alpha}
\donttest{result$wpost}

# With prior on structure among predictors
data(linearrelation)
\donttest{result<-icmm(Y=Y, X=X, b.start=initbeta, family="gaussian", 
        ising.prior=TRUE, structure=linearrelation, 
        estalpha=FALSE, alpha=0.5, maxiter=100)}
\donttest{result$coef}
\donttest{result$iterations}
\donttest{result$alpha}
\donttest{result$wpost}

# Binary logistic regression model
data(simBinomial)
Y<-as.matrix(simBinomial[,1])
X<-as.matrix(simBinomial[,-1])
p<-dim(X)[2]
# Obtain initial values from lasso
data(initbetaBinomial)
initbeta<-as.matrix(initbetaBinomial)
\donttest{result<-icmm(Y=Y, X=X, b0.start=0, b.start=initbeta, family="binomial", 
        ising.prior=TRUE, structure=linearrelation, estalpha=FALSE, 
        alpha=0.5, maxiter=100)}
\donttest{result$coef}
\donttest{result$iterations}
\donttest{result$alpha}
\donttest{result$wpost}

# Cox's model
data(simCox)
Y<-as.matrix(simCox[,1])
event<-as.matrix(simCox[,2])
X<-as.matrix(simCox[,-(1:2)])
# Obtain initial values from lasso
data(initbetaCox)
initbeta<-as.matrix(initbetaCox)
\donttest{result <- icmm(Y=Y, X=X, event=event, b.start=initbeta, family="cox",
          ising.prior=TRUE, structure=linearrelation, estalpha=FALSE,
          alpha=0.5, maxiter=100)}
\donttest{result$coef}
\donttest{result$iterations}
\donttest{result$alpha}
\donttest{result$wpost}
}
