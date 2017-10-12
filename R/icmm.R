icmm<-function(Y, X, event, b0.start, b.start, family="gaussian", ising.prior=FALSE, structure, estalpha=FALSE, alpha=0.5, maxiter=100)
{
  # check function arguments and inputs
  if (missing(b0.start))
  {
    b0.start<-0
  }
  if (family!="gaussian" & family!="binomial" & family!="cox")
  {
    stop("family must be 'gaussian', 'binomial', or 'cox'.")
  }
  if (family=="cox")
  {
    if (missing(event))
    {
      stop("event argument is missing.")
    }
  }
  if (ising.prior==TRUE)
  {
    if (missing(structure))
    {
      stop("structure argument is missing.")
    }
  }
  n<-dim(X)[1]
  p<-dim(X)[2]
  nY<-dim(Y)[1]
  if (nY!=n)
  {
    stop("Y and X must have same number of samples.")
  }
  pb.start<-dim(b.start)[1]
  if (pb.start!=p)
  {
    stop("b.start must have same number of predictos as in X.")
  }

  # initiate all variables
  niter<-0
  tol<-10^(-6)
  epsilon<-10^(-50)
  diff<-10
  alpha<-alpha

  if (family=="gaussian")
  {
    # standardize X and centering Y
    X<-scale(X,center=T,scale=T)
    x.mean<-attr(X, 'scaled:center')
    x.sd<-attr(X, 'scaled:scale')
    Y<-scale(Y,center=T,scale=F)
    y.mean<-attr(Y,'scaled:center')
    b<-b.start
    b<-b*x.sd
    sigma<-get.sigma(Y,X,b,alpha)
  }
  if (family=="binomial")
  {
    b0<-b0.start
    b<-b.start
  }
  if (family=="cox")
  {
    b<-b.start
    time<-sort(unique(Y))
    ntime<-length(time)
    # sum of event_i where y_i =time_k
    sumevent<-rep(0, ntime)
    for(j in 1:ntime)
    {
      sumevent[j]<-sum(event[Y[,1]==time[j]])
    }
  }

  if (ising.prior==TRUE)
  {
    edgeind<-sort(unique(structure[,1]))
  }
  else
  {
    w<-0.5
  }

  # The iterative procedure starts here
  while (niter<maxiter && diff>tol)
  {
    b_old<-b
    nzind_old<-which(b_old!=0)

    if (family=="gaussian")
    {
      if (ising.prior==TRUE)
      {
        # Updata a and b (hyperpara[1]=a, hyperparam[2]=b)
        hyperparam<-get.ab(b, structure, edgeind)
      }
      # Update beta
      for (j in 1:p)
      {
        Yres<-Y-X%*%b+X[,j]*b[j,1]
        sxy<-t(Yres)%*%X[,j]
        ssx<-sum(X[,j]^2)
        SS<-sqrt(n-1)*sxy/(sigma*ssx)
        if (ising.prior==FALSE)
        {
          b[j,1]<-get.beta(SS, w, alpha, sigma/sqrt(n-1))
        }
        else
        {
          wpost<-get.wpost(SS, b, alpha, hyperparam, structure, edgeind, j)
          b[j,1]<-get.beta.ising(SS, wpost, alpha, sigma/sqrt(n-1))
        }
      }
      # Update sigma
      sigma<-get.sigma(Y,X,b,alpha)
      # Update alpha
      if (estalpha==TRUE)
      {
        alpha<-get.alpha(b, 1/(sqrt(n-1)*sum(abs(b))/sigma))
      }
      # Update w
      if (ising.prior==FALSE)
      {
        w<-get.wprior(b)
      }
      diff<-sum((b-b_old)^2)/ifelse(sum(b_old^2)==0, epsilon, sum(b_old^2))
    }

    if (family!="gaussian")  # for binomial and cox model
    {
      # get pseudodata
      if (family=="binomial")
      {
        pseudodata<-get.pseudodata.binomial(Y, X, b0, b, niter)
      }
      if (family=="cox")
      {
        pseudodata<-get.pseudodata.cox(Y, X, event, b, time, ntime, sumevent)
        b0<-0
      }
      z<-pseudodata$z
      sigma<-matrix(sqrt(pseudodata$sigma2), ncol=1)

      if (ising.prior==TRUE)
      {
        # Updata a and b (hyperpara[1]=a, hyperparam[2]=b)
        hyperparam<-get.ab(b, structure, edgeind)
      }
      # Update beta
      X2<-X^2
      X2S2<-t(1/(sigma^2))%*%X2
      for (j in 1:p)
      {
        Zres<-z-b0-X%*%b+X[,j]*b[j,1]
        ZS2<-Zres/(sigma^2)
        XZS2<-t(ZS2)%*%X[,j]
        SS<-XZS2/sqrt(X2S2[1,j])
        if (ising.prior==FALSE)
        {
          b[j,1]<-get.beta(SS, w, alpha, 1/sqrt(X2S2[1,j]))
        }
        else
        {
          wpost<-get.wpost(SS, b, alpha, hyperparam, structure, edgeind, j)
          b[j,1]<-get.beta.ising(SS, wpost, alpha, 1/sqrt(X2S2[1,j]))
        }
      }

      # get intercept term (b0) for binomial case
      if (family=="binomial")
      {
        weight<-1/(sigma^2)
        tmp<-matrix(apply(X, 2,weighted.mean,w=weight), nrow=1)
        b0<-as.numeric(weighted.mean(z, weight)-tmp%*%b)
      }

      # Update alpha
      if (estalpha==TRUE)
      {
        alpha<-get.alpha(b, 1/(sqrt(X2S2)%*%abs(b)))
      }
      # Update w
      if (ising.prior==FALSE)
      {
        w<-get.wprior(b)
      }

      # check for activeset convergence
      nzind<-which(b!=0)
      if (setequal(nzind_old, nzind))
      {
        diff<-0
      }
    }

    niter<-niter+1
    if (niter%%10==0)
    {
      print(paste("iterations: ", niter, ";  relative difference: ", diff))
    }
  }

  # Calculate zeta (local posterior probability)
  zeta<-rep(0,p)
  if (family=="gaussian")
  {
    for (j in 1:p)
    {
      Yres<-Y-X%*%b+X[,j]*b[j,1]
      sxy<-t(Yres)%*%X[,j]
      ssx<-sum(X[,j]^2)
      SS<-sqrt(n-1)*sxy/(sigma*ssx)

      if (ising.prior==FALSE)
      {
        zeta[j]<-get.zeta(SS, w, alpha)
      }
      else
      {
        zeta[j]<-get.zeta.ising(SS, b, alpha, hyperparam, structure, edgeind,j)
      }
    }
  }

  if (family!="gaussian")
  {
    nz<-sum(b!=0)
    X2<-X^2
    X2S2<-t(1/(sigma^2))%*%X2
    for (j in 1:p)
    {
      Zres<-z-b0-X%*%b+X[,j]*b[j,1]
      ZS2<-Zres/(sigma^2)
      XZS2<-t(ZS2)%*%X[,j]
      SS<-XZS2/sqrt(X2S2[1,j])

      if (ising.prior==FALSE)
      {
        zeta[j]<-get.zeta(SS, w, alpha)
      }
      else
      {
        zeta[j]<-get.zeta.ising(SS, b, alpha, hyperparam, structure, edgeind,j)
      }
    }
  }

  if (estalpha==TRUE)
  {
    print(paste("Optimum alpha: ", alpha, sep=""))
  }
  print(paste("Number of iterations: ", niter, sep=""))

  # transform back into the original model
  if (family=="gaussian")
  {
    output<-matrix(rep(0,p+1), ncol=1)
    # intercept
    output[1,1]<-y.mean-sum(x.mean*b/x.sd)
    output[2:(p+1),1]<-b/x.sd
  }

  if (family=="binomial")
  {
    output<-matrix(rep(0,p+1), ncol=1)
    # intercept
    output[1,1]<-b0
    output[2:(p+1),1]<-b
  }

  if (family=="cox")
  {
    output<-matrix(b[,1], ncol=1)
  }

  result<-list(as.vector(output), niter, alpha, zeta)
  names(result) <- c("coef", "iterations", "alpha", "postprob")
  return(result)
}



