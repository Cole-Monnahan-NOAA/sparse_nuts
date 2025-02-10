# remotes::install_github("stan-dev/posteriordb-r")
#library(tidyverse)
# library(posteriordb)
# my_pdb <- pdb_github()
# po <- 'gp_pois_regr'
# dat <- pdb_data(po)
# code <- stan_code(po)
# ref <- reference_posterior_draws(po)
# df <- lapply(ref[1:4], \(x) as.data.frame(x)) |> bind_rows()
# saveRDS(dat, file='models/gp_pois_regr/dat.RDS')
# saveRDS(df, file='models/gp_pois_regr/posterior.RDS')
# writeLines(code, con ='models/gp_pois_regr/code.stan')
# dat <- readRDS('models/gp_pois_regr/dat.RDS')
# dat <- list(x=dat$x, k=dat$k)

# reference draws from posteriordb
df <- readRDS('models/gp_pois_regr/ref_samples.RDS')
cor <- cor(df)
corrplot::corrplot(cor)
#pairs(df, upper.panel = NULL, pch='.')


library(RTMB)
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
  return(-lp) # note the negative log posterior
}
obj0 <- MakeADFun(func=func, parameters=pars, random=NULL, silent=TRUE)
stopifnot(obj0$fn()==func(pars))

# will escape warmup with seed=1, will not with seed=2
library(StanEstimators)
test <- stan_sample(fn=function(x) -obj0$fn(x),
                    grad_fun = function(x) -obj0$gr(x), num_chains=1, seed = 1,
                    num_samples=1000, num_warmup=1000, par_inits=obj0$par)

obj <- MakeADFun(func=func, parameters=pars, random='f_tilde', silent=TRUE)
opt <- with(obj, nlminb(par, fn, gr))
exp(opt$par)


# run rstan to get f_tilde and check that it matches the RTMB model
library(rstan)
dat$N <- 11
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanfit <- stan('models/gp_pois_regr/gp_pois_regr.stan', data=dat,  chains=4, cores=4)
df2 <- as.data.frame(stanfit)
## diff two lps to ditch the constants
par2 <- df2[1,] |> as.numeric()
pars1 <- list(logrho=log(par2[1]), logalpha=log(par2[2]),
              f_tilde=par2[3:13])
par2 <- df2[20,] |> as.numeric()
pars2 <- list(logrho=log(par2[1]), logalpha=log(par2[2]),
              f_tilde=par2[3:13])
## matches nicely
abs(-(func(pars1)-func(pars2)) - (df2$lp__[1]-df2$lp__[20]))


## double check that pars are sampled the same
library(adnuts)
fit <- sample_sparse_tmb(obj, iter=5000, warmup=1000, chains=1,
                         cores=1, globals=list(dat=dat), control=list(adapt_delta=.95),
                         seed = 2)

pairs_admb(fit, pars=1:6, order='slow')
pairs_admb(fit, pars=1:3)
plot_uncertainties(fit)
plot_sampler_params(fit)

# check that pars are the same, first have to massage TMB output to match Stan's
post <- as.data.frame(fit)
# for some reason they report f and not f_tilde?
frep <- t(apply(post, 1, \(x) obj$report(x)$f))
df.tmb <- data.frame(rho=exp(post$logrho), alpha=exp(post$logalpha), frep)
names(df.tmb)=names(df)

library(tidyr); library(dplyr); library(ggplot2)
df.all <- bind_rows(
  pivot_longer(df, cols=everything()) |> mutate(platform='posteriordb'),
  pivot_longer(df.tmb, cols=everything()) |> mutate(platform='TMB'))
ggplot(df.all, aes(x=platform, y=value)) + facet_wrap('name') + geom_violin()

