library(rstan)

setwd(here::here('models/8schools'))
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
obj <- MakeADFun(f, pars, random=c("eta"), silent=TRUE)

## run longer chains
library(SparseNUTS)
library(StanEstimators)
fit <- sample_snuts(obj, iter=2000, warmup=200, chains=4,
                         cores=4, globals=list(schools_dat=schools_dat))

## compare correlations and marginal sds
post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))

obj$par |> length()
obj$env$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

minsd <- apply(post, 2, sd) |> min()
maxsd <- apply(post, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd



pairs(fit, pars=1:6, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)

setwd(here::here())

