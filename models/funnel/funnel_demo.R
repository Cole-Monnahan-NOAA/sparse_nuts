# Funnel example ported to RTMB from
# https://mc-stan.org/docs/cmdstan-guide/diagnose_utility.html#running-the-diagnose-command

#install.packages('RTMB')
library(RTMB)
pars <- list(y=-1.12, x=rep(0,9)) # initialize and declare dimensions
## the (negative) posterior density as a function in R
f <- function(pars){
  getAll(pars)
  lp <- dnorm(y, 0, 3, log=TRUE) + # prior
    sum(dnorm(x, 0, exp(y/2), log=TRUE)) # likelihood
  return(-lp) # TMB expects negative log posterior
}
# this is the TMB object with autodiff lp and gradients of the marginal
obj <- RTMB::MakeADFun(f, pars, random='x', silent=TRUE)
obj$fn()
obj$gr()

### Now SNUTS
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
library(adnuts)
library(StanEstimators)
fit <- sample_sparse_tmb(obj, iter=1000, warmup=200, seed=1213)
fit$mle$Q    # joint precision
fit$mle$Qinv # joint covariance
post <- as.data.frame(fit)
# hasn't recovered the prior b/c it's not converged, particularly
# for small y values
hist(post$y, freq=FALSE)
lines(x<-seq(-10,10, len=200), dnorm(x,0,3))
abline(v=fit$mle$est[1], col=2, lwd=2)

# do ELA so there is only one parameter: y = 2*log SD
fit.ela <- sample_sparse_tmb(obj, iter=2000, warmup=200, laplace=TRUE,
                             seed=12312)
# you just get the prior back b/c the Laplace approximation is
# accurate
post.ela <- as.data.frame(fit.ela)
hist(post.ela$y, freq=FALSE, breaks=30)
lines(x<-seq(-10,10, len=200), dnorm(x,0,3))
