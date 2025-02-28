# source('radon.R')
#  dat <- list(y=y, x=x, u=u, N=N, J=J, holdout=holdout, county=county)
#  saveRDS(dat, file='models/radon/dat.RDS')

setwd(here::here('models/radon'))

dat <- readRDS('dat.RDS')

library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)

fit.stan <- stan(file='radon.stan', data = dat)
library(shinystan)
launch_shinystan(fit.stan)

library(RTMB)
dat.radon <- dat
pars <- list(a=rep(0, dat$J), mu_a=1,  logsigma_a=1, logsigma_y=1)
f <- function(pars){
  getAll(pars,dat)
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
f(pars)
obj <- MakeADFun(f, pars, random='a', silent=TRUE)
obj$par |> length()
obj$env$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))


post <- as.data.frame(fit.stan)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))
p1 <- post[1,-390]
p2 <- post[1000,-390]
diff.stan <- post[1,390]-post[1000,390]
p1$sigma_a <- log(p1$sigma_a)
p1$sigma_y <- log(p1$sigma_y)
p2$sigma_a <- log(p2$sigma_a)
p2$sigma_y <- log(p2$sigma_y)

diff.tmb <- obj$report(as.numeric(p1))$lp-
  obj$report(as.numeric(p2))$lp

diff.stan - diff.tmb

library(adnuts)
fit <- sample_sparse_tmb(obj, iter=2000, warmup=1000, chains=4,
                         cores=4, globals=list(dat=dat), metric='unit')

pairs_admb(fit, pars=1:7, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)


post.tmb <- as.data.frame(fit)
minsd <- apply(post.tmb, 2, sd) |> min()
maxsd <- apply(post.tmb, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd
