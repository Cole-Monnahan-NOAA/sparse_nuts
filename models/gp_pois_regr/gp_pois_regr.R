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
setwd(here::here('models/gp_pois_regr'))

# reference draws from posteriordb
df <- readRDS('ref_samples.RDS')
cor <- cor(df)
corrplot::corrplot(cor)
#pairs(df, upper.panel = NULL, pch='.')
dat <- list(x = c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10),
            k = c(40, 37, 29, 12, 4, 3, 9, 19, 77, 82, 33))
pars <- list(logrho=1, logalpha=1,  f_tilde=rep(0,11))

library(TMB)
compile("gp_pois_regr.cpp")
dyn.unload('gp_pois_regr')
dyn.load('gp_pois_regr')
obj <- MakeADFun(data=dat, parameters=pars,  DLL='gp_pois_regr')
#obj$report()$lp -lp

# run rstan to get f_tilde and check that it matches the RTMB model
library(rstan)
dat$N <- 11
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanfit <- stan('gp_pois_regr.stan', data=dat,  chains=4, cores=4)
df2 <- as.data.frame(stanfit)
## diff two lps to ditch the constants
par2 <- df2[1,] |> as.numeric()
pars1 <- list(logrho=log(par2[1]), logalpha=log(par2[2]),
              f_tilde=par2[3:13])
par2 <- df2[20,] |> as.numeric()
pars2 <- list(logrho=log(par2[1]), logalpha=log(par2[2]),
              f_tilde=par2[3:13])
## matches nicely
abs(-(obj$fn(unlist(pars1))-obj$fn(unlist(pars2))) - (df2$lp__[1]-df2$lp__[20]))
#abs(-(func(pars1)-func(pars2)) - (df2$lp__[1]-df2$lp__[20]))

## longer chains in TMB
obj <- MakeADFun(data=dat, parameters=pars,  DLL='gp_pois_regr', random='f_tilde')
setwd(here::here('models/gp_pois_regr'))
saveRDS(obj, 'obj.gp_pois_regr.RDS')
setwd(here::here())
fit <- sample_sparse_tmb(obj, iter=8000, warmup=200, chains=4,
                         thin=4,
                         seed=1, control=list(adapt_delta=.99),
                         cores=4, globals=list(dat=dat))
pairs_admb(fit, pars=1:5, order='slow')

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

### -------------------------
# Old RTMB code. This works but is unstable due to the Chol call.
# No graceful way to catch this so StanEstimators crashes. So
# switched to TMB above.

# library(RTMB)
# func <- function(pars){
#   RTMB::getAll(pars,dat)
#   alpha <- exp(logalpha)
#   rho <- exp(logrho)
#   cov <- matrix(NA, 11,11)
#   for(i in 1:11){
#     for(j in 1:11){
#       # https://mc-stan.org/docs/functions-reference/matrix_operations.html#exponentiated-quadratic-kernel
#       cov[i,j] <- alpha^2*exp(-(x[i]-x[j])^2/(2*rho^2))
#     }
#   }
#   cov <- cov+diag(1e-10,11)
#   L_cov <- t(chol(cov)) # upper tri chol
#   f <- as.numeric(L_cov %*% f_tilde) # log-predicted lambda
#   lp <-
#     # priors
#     dgamma(rho, shape=25,scale=1/4, log=TRUE)+
#     dnorm(alpha, 0,2,log=TRUE)+
#     sum(dnorm(f_tilde,0,1,log=TRUE)) + # hyperdistribution
#     logrho + logalpha + # Jacobians
#     sum(dpois(k, exp(f), log=TRUE)) # likelihood
#   REPORT(f)
#   return(-lp) # note the negative log posterior
# }

# #stopifnot(obj0$fn()==func(pars))
# obj <- MakeADFun(func=func, parameters=pars, random='f_tilde', silent=TRUE)
# opt <- with(obj, nlminb(par, fn, gr))


# # will escape warmup with seed=1, will not with seed=2
# library(StanEstimators)
# obj0 <- MakeADFun(func=func, parameters=pars, random=NULL, silent=TRUE)
# test <- stan_sample(fn=function(x) -obj0$fn(x),
#                     grad_fun = function(x) -obj0$gr(x), num_chains=1, seed = 1,
#                     num_samples=1000, num_warmup=1000, par_inits=obj0$par)

