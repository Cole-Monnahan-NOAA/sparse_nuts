
if('TMB' %in% .packages()) detach(package:TMB)

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
obj.kilpisjarvi <- MakeADFun(f, pars, random=NULL, silent=TRUE)


radon_dat <- readRDS('models/radon/dat.RDS')
pars <- list(a=rep(0, radon_dat$J), mu_a=1,  logsigma_a=1, logsigma_y=1)
f <- function(pars){
  getAll(pars,radon_dat)
  y_hat <- a[county]
  sigma_a <- exp(logsigma_a)
  sigma_y <- exp(logsigma_y)
  lp <-
    dnorm(mu_a, 0,1, log=TRUE) +
    sum(dnorm(a, mu_a, sigma_a, log=TRUE)) +
    sum(dnorm(y, y_hat, sigma_y, log=TRUE)) +
    logsigma_a + logsigma_y
  REPORT(lp)
  return(-lp)
}
obj.radon <- MakeADFun(f, pars, random='a', silent=TRUE, DLL='radon')


library(dsem)
inputs <- readRDS('models/causal/inputs.RDS')
obj.causal <- with(inputs, dsemRTMB( sem = sem,
                    tsdata = tsdata,
                    family = family,
                    log_prior = log_prior,
                    control = dsem_control(run_model=FALSE, use_REML = FALSE) ))$obj

