setwd(here::here('models/accel_gp'))

library(MASS)
library(brms)
library(rstan)

code = make_stancode(formula = bf(accel ~ gp(times, k=40, c=3/2), sigma ~ gp(times, k=20, c=3/2)), data=mcycle)
data = make_standata(formula = bf(accel ~ gp(times, k=40, c=3/2), sigma ~ gp(times, k=20, c=3/2)), data=mcycle)

write(code, "accel_gp.stan")
#stan_rdump(names(data), "accel_gp.dat", envir = list2env(data))
saveRDS(data, file='dat.RDS')

dat <- readRDS('dat.RDS')


library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

fit.stan <- stan(file='accel_gp.stan', data = dat)
library(shinystan)
launch_shinystan(fit.stan)

pars <- list(alpha=1, beta=1, logsigma=0)

library(RTMB)

f <- function(pars){
  getAll(pars,dat)
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
obj <- MakeADFun(f, pars, random=NULL, silent=TRUE)

opt <- with(obj, nlminb(par,fn,gr))

post <- as.data.frame(fit.stan)
cors <- cor(post)
corrplot::corrplot(cors, type='lower')
p1 <- post[1,-4]
p2 <- post[1000,-4]
diff.stan <- post[1,4]-post[1000,4]
p1$sigma <- log(p1$sigma)
p2$sigma <- log(p2$sigma)

diff.tmb <- obj$report(as.numeric(p1))$lp-
  obj$report(as.numeric(p2))$lp

diff.stan - diff.tmb

library(SparseNUTS)
fit <- sample_snuts(obj, iter=2000, warmup=1000, chains=4,
                         cores=4, globals=list(dat=dat))

pairs(fit, pars=1:7, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)


post.tmb <- as.data.frame(fit)
apply(post.tmb, 2, sd)


# N <- 62
# x <-
#   c(3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3968,
#     3969, 3970, 3971, 3972, 3973, 3974, 3975, 3976, 3977, 3978, 3979, 3980, 3981, 3982, 3983, 3984, 3985, 3986,
#     3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997, 3998, 3999, 4000, 4001, 4002, 4003, 4004,
#     4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012, 4013)
# xpred <- 2016
# y <-
#   c(8.3, 10.9, 9.4, 8.1, 8.1, 7.7, 8.6, 9.1, 11, 10.1, 7.6, 8.8, 8.3, 7.2, 9.3, 8.8, 7.6, 10.5, 11, 8.9,
#     11.3, 10, 10.1, 6.4, 8.2, 8.4, 9.5, 9.9, 10.6, 7.6, 7.7, 8.1, 8.4, 9.7, 9.5, 7.3, 10.3, 9.6, 10.3, 9.8, 9, 9.1,
#     9.5, 8.7, 9.9, 10.5, 9.4, 9, 9, 9.7, 11.4, 10.7, 10.1, 10.8, 10.4, 10.3, 8.8, 9.8, 8.8, 10.8, 8.6, 11.1)
# pmualpha <- 9.31290322580645
# psalpha <- 100
# pmubeta <- 0
# psbeta <- 0.0333333333333333
# dat <- list(N=N, x=x, xpred=xpred, y=y, pmualpha=pmualpha, psalpha=psalpha, pmubeta=pmubeta, psbeta=psbeta)
# saveRDS(dat, 'dat.RDS')

dat <- readRDS('dat.RDS')


library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

fit.stan <- stan(file='accel_gp.stan', data = dat)
library(shinystan)
launch_shinystan(fit.stan)

pars <- list(alpha=1, beta=1, logsigma=0)

library(RTMB)

f <- function(pars){
  getAll(pars,dat)
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
obj <- MakeADFun(f, pars, random=NULL, silent=TRUE)

opt <- with(obj, nlminb(par,fn,gr))

post <- as.data.frame(fit.stan)
cors <- cor(post)
corrplot::corrplot(cors, type='lower')
p1 <- post[1,-4]
p2 <- post[1000,-4]
diff.stan <- post[1,4]-post[1000,4]
p1$sigma <- log(p1$sigma)
p2$sigma <- log(p2$sigma)

diff.tmb <- obj$report(as.numeric(p1))$lp-
  obj$report(as.numeric(p2))$lp

diff.stan - diff.tmb

library(SparseNUTS)
fit <- sample_snuts(obj, iter=2000, warmup=1000, chains=4,
                         cores=4, globals=list(dat=dat))

pairs(fit, pars=1:7, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)


post.tmb <- as.data.frame(fit)
apply(post.tmb, 2, sd)
