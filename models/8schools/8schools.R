library(rstan)


schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = 'schools.stan', data = schools_dat)


library(RTMB)

pars <- list(mu=0, logtau=0, eta=rep(1,8))

f <- function(pars){
  getAll(schools_dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # prior
         sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
         logtau                          # jacobian
  return(-lp)
}
obj.schools <- MakeADFun(f, pars, random=c("eta"))
opt <- TMBhelper::fit_tmb(obj, getJointPrecision=TRUE)

library(adnuts)
library(StanEstimators)
fit <- sample_sparse_tmb(obj.schools, iter=2000, warmup=100, cores=1, chains=4)
pairs_admb(fit)
