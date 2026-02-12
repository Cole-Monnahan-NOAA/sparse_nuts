library(RTMB)

setwd(here::here('models/funnel'))

library(rstan)
fit.stan <- stan('funnel.stan', model_name='funnel',
                 data=list(), cores=1)



pars <- list(y=0, x=rep(1,9))

library(RTMB)
## define the (negative) posterior density as a function in R
f <- function(pars){
  getAll(pars)
  lp <- dnorm(y, 0, 3, log=TRUE) + # prior
    sum(dnorm(x, 0, exp(y/2), log=TRUE)) # likelihood
  REPORT(lp)
  return(-lp) # TMB expects negative log posterior
}
obj <- RTMB::MakeADFun(f, pars, random='x', silent=TRUE)
obj$fn()


# check TMB matches Stan
post <- as.data.frame(fit.stan)
p1 <- post[1,-11]
p2 <- post[1000,-11]
diff.stan <- post$lp__[1]-post$lp[1000]
diff.tmb <- obj$report(as.numeric(p1))$lp-
  obj$report(as.numeric(p2))$lp
diff.stan - diff.tmb # zero if the model matches

## run longer chains
library(SparseNUTS)
library(StanEstimators)
fit <- sample_snuts(obj, iter=2000, warmup=200, chains=4,
                         cores=1)
pairs(fit, pars=1:4)
post.tmb <- as.data.frame(fit)

# do ELA
fit.ela <- sample_sparse_tmb(obj, iter=2000, warmup=200, chains=4,
                             cores=1, laplace=TRUE)
pairs(fit.ela)

# you just get the prior back b/c the Laplace approximation is
# accurate
post.ela <- extract_samples(fit.ela)
hist(post.ela$y, freq=FALSE, breaks=30)
lines(x<-seq(-10,10, len=200), dnorm(x,0,3))
