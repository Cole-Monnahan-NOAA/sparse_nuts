
library(RTMB)
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(mu=0, logtau=0, eta=rep(1,8))

f <- function(pars){
  RTMB::getAll(schools_dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # prior
    sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
    logtau                          # jacobian
  return(-lp)
}
obj.schools <- RTMB::MakeADFun(func=f, parameters=pars,
                               random="eta", silent=TRUE)



diamonds_dat <- readRDS(file='models/diamonds/dat.RDS')
Kc <- diamonds_dat$K - 1;
Xc <- matrix(NA, nrow=diamonds_dat$N, ncol=Kc) # centered version of X without an intercept
means_X <- rep(NA, Kc) # column means of X before centering
for(i in 2 : diamonds_dat$K) {
  means_X[i - 1] = mean(diamonds_dat$X[,i]);
  Xc[,i - 1] = diamonds_dat$X[,i] - means_X[i - 1];
}
diamonds_dat$means_X <- means_X
diamonds_dat$Xc <- Xc
pars <- list(b=rep(0, Kc),## population-level effects
             Intercept=1,   ## temporary intercept for centered predictors
             logsigma=0) ## residual SD

## define the (negative) posterior density as a function in R
f <- function(pars){
  getAll(pars,diamonds_dat)
  ## transformed parameters
  sigma <- exp(logsigma)
  lp <- logsigma + ## Jacobian
    ## priors
    sum(dnorm(b, 0, 1, log=TRUE)) +
    sum(dt((Intercept-8)/10, 3,log=TRUE) )+ # dropped constants
    sum(dt(sigma/10,df=3, log=TRUE))
  ## likelihood including all constants
  mu_hat <- Intercept + as.numeric(Xc%*%b)
  if (!prior_only)  lp <- lp+sum(dnorm(Y, mean=mu_hat, sd=sigma, log=TRUE))
  ## generated quantities
  REPORT(mu_hat)
  b_Intercept <- Intercept - sum(means_X*b)
  REPORT(b_Intercept)
  REPORT(lp)
  return(-lp) # TMB expects negative log posterior
}
obj.diamonds <- RTMB::MakeADFun(f, pars, random='b', silent=TRUE)


dat <- list(x = c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10),
            k = c(40, 37, 29, 12, 4, 3, 9, 19, 77, 82, 33))
pars <- list(logrho=1, logalpha=1,  f_tilde=rep(0,11))
func <- function(pars){
  getAll(pars,dat)
  alpha <- exp(logalpha)
  rho <- exp(logrho)
  cov <- matrix(NA, 11,11)
  for(i in 1:11){
    for(j in 1:11){
      # https://mc-stan.org/docs/functions-reference/matrix_operations.html#exponentiated-quadratic-kernel
      cov[i,j] <- alpha^2*exp(-(x[i]-x[j])^2/(2*rho^2))
    }
  }
  cov <- cov+diag(1e-10,11)
  L_cov <- t(chol(cov)) # upper tri chol
  f <- as.numeric(L_cov %*% f_tilde) # log-predicted lambda
  lp <-
    # priors
    dgamma(rho, shape=25,scale=1/4, log=TRUE)+
    dnorm(alpha, 0,2,log=TRUE)+
    sum(dnorm(f_tilde,0,1,log=TRUE)) + # hyperdistribution
    logrho + logalpha + # Jacobians
    sum(dpois(k, exp(f), log=TRUE)) # likelihood
  REPORT(f)
  return(-lp)
}

obj.gp_pois <- MakeADFun(func=func, parameters=pars, random='f_tilde', silent=TRUE)




kilpisjarvi_dat <- readRDS('models/kilpisjarvi/dat.RDS')
pars <- list(alpha=1, beta=1, logsigma=0)
f <- function(pars){
  getAll(pars,kilpisjarvi_dat)
  sigma <- exp(logsigma)
  lp <-
    dnorm(alpha, pmualpha, psalpha, log=TRUE) +
    dnorm(beta, pmubeta, psbeta, log=TRUE) +
    sum(dnorm(y, alpha+beta*x, sigma, log=TRUE)) +
    logsigma
  REPORT(lp)
  return(-lp)
}

f(pars)
obj.kilpisjarvi <- MakeADFun(f, pars, random=NULL, silent=TRUE)
