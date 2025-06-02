setwd(here())

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



irt_2pl_dat <- readRDS('models/irt_2pl/dat.RDS')
pars <- list(logsigma_theta=1, theta=rep(1, irt_2pl_dat$J),
             logsigma_a=-0.8, loga=rep(.5, irt_2pl_dat$I),
             mu_b=1, logsigma_b=1, b=rep(.5, irt_2pl_dat$I))
func <- function(pars){
  getAll(irt_2pl_dat, pars)
  sigma_theta <- exp(logsigma_theta)
  sigma_a <- exp(logsigma_a)
  sigma_b <- exp(logsigma_b)
  a <- exp(loga)
  lp <-
    #  jacobian
    logsigma_a+logsigma_b+  logsigma_theta+
    # priors
    dcauchy(sigma_theta, 0, 2, log=TRUE)+
    sum(dnorm(theta, 0, sigma_theta, log=TRUE))+
    dcauchy(sigma_a, 0, 2, log=TRUE) +
    sum(dnorm(loga, 0, sigma_a, log=TRUE))+
    dnorm(mu_b, mean = 0, sd=5, log=TRUE) +
    dcauchy(sigma_b, 0, 2, log=TRUE)+
    # hyperdistribution
    sum(dnorm(b, mu_b, sigma_b, log=TRUE))
  # likelihood
  for(i in 1:I) {
    pred <- plogis(a[i]*(theta-b[i]))
    lp <- lp + sum(RTMB:::Term(dbinom(y[i,], size=1, prob=pred, log=TRUE)))
  }
  REPORT(lp)
  return(-lp)
}
obj.irt_2pl <- MakeADFun(func, pars, random=c('theta', 'loga','b'))

func_nc <- function(pars){
  getAll(irt_2pl_dat, pars)
  sigma_theta <- exp(logsigma_theta)
  sigma_a <- exp(logsigma_a)
  sigma_b <- exp(logsigma_b)
  a <- exp(sigma_a*loga)
  lp <-
    #  jacobian
    logsigma_a+logsigma_b+  logsigma_theta+
    # priors
    dcauchy(sigma_theta, 0, 2, log=TRUE)+
    sum(dnorm(theta, 0, sigma_theta, log=TRUE))+
    dcauchy(sigma_a, 0, 2, log=TRUE) +
    sum(dnorm(loga, 0, 1, log=TRUE))+
    dnorm(mu_b, mean = 0, sd=5, log=TRUE) +
    dcauchy(sigma_b, 0, 2, log=TRUE)+
    # hyperdistribution
    sum(dnorm(b, mu_b, sigma_b, log=TRUE))
  # likelihood
  for(i in 1:I) {
    pred <- plogis(a[i]*(theta-b[i]))
    lp <- lp + sum(RTMB:::Term(dbinom(y[i,], size=1, prob=pred, log=TRUE)))
  }
  REPORT(lp)
  return(-lp)
}
obj.irt_2pl_nc <- MakeADFun(func_nc, pars, random=c('theta', 'loga','b'))




# local copy of ecostate::rk4sys
rk4sys <- function (f, a, b, y0, n, Pars, ...){
  m <- length(y0)
  h <- (b - a)/n
  x <- seq(a + h, b, by = h)
  y <- matrix(0, nrow = n, ncol = m)
  k1 <- h * f(a, y0, Pars, ...)
  k2 <- h * f(a + h/2, y0 + k1/2, Pars, ...)
  k3 <- h * f(a + h/2, y0 + k2/2, Pars, ...)
  k4 <- h * f(a + h, y0 + k3, Pars, ...)
  y[1, ] <- y0 + k1/6 + k2/3 + k3/3 + k4/6
  for (i in seq_len(n - 1)) {
    k1 <- h * f(x[i], y[i, ], Pars, ...)
    k2 <- h * f(x[i] + h/2, y[i, ] + k1/2, Pars, ...)
    k3 <- h * f(x[i] + h/2, y[i, ] + k2/2, Pars, ...)
    k4 <- h * f(x[i] + h, y[i, ] + k3, Pars, ...)
    y[i + 1, ] <- y[i, ] + k1/6 + k2/3 + k3/3 + k4/6
  }
  x <- c(a, x)
  y <- rbind(y0, y)
  return(list(x = x, y = y))
}
f <- identity
dzdt <- function(Time, State, Pars){
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  d1_dt = f(Pars$ln_alpha) * State[1] - f(Pars$ln_beta) * State[1] * State[2]
  d2_dt = -1 * f(Pars$ln_gamma) * State[2] + f(Pars$ln_delta) * State[1] * State[2]
  return( c(d1_dt, d2_dt) )
}

get_nll <- function(pars){
  getAll(pars)
  zhat_ti <- array(0, dim = dim(dat) )
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  # Initialize
  zhat_ti[1,] <- ln_z0
  for( t_index in 2:nrow(dat) ){
    proj <- rk4sys(f=dzdt, a=0, b=1, n=10,
                   Pars=pars, y0=z_ti[t_index-1,])
    zhat_ti[t_index,] <- proj$y[nrow(proj$y),]
  }
  lp <- ln_sigma + ln_tau + # jacobians
    sum(dnorm(dat, mean=z_ti, sd=exp(ln_sigma), log=TRUE)) +
    sum(dnorm(z_ti, mean=zhat_ti, sd=exp(ln_tau), log=TRUE)) +
    # informative prior
    dgamma(x=exp(ln_sigma), scale=.03, shape=9.5, log=TRUE)
  ADREPORT( zhat_ti )
  REPORT(lp)
  return(-lp)
}
# https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html
hares <- c(30, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22, 25.4,
           27.1, 40.3, 57, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7)
lynx <- c(4, 6.1, 9.8, 35.2, 59.4, 41.7, 19, 13, 8.3, 9.1, 7.4,
          8, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6)
dat <- cbind(X=hares/10, Y=lynx/10)
pars <- list(ln_alpha=1, ln_beta=1, ln_gamma=1,
             ln_delta=1, ln_z0=c(0.1,0.1),
             ln_sigma=1, ln_tau=1,
             z_ti = dat) #ifelse( is.na(dat), rep(1,nrow(dat))%o%colMeans(dat,na.rm=TRUE), dat
globals.lynx <- list(dat=dat, rk4sys=rk4sys, dzdt=dzdt, f=f)
get_nll(pars)
obj.lynx <- MakeADFun( func = get_nll, silent=TRUE,
                  parameters = pars,
                  random = "z_ti" )
