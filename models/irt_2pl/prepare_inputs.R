setwd(here::here('models/irt_2pl'))
#
# # remotes::install_github("stan-dev/posteriordb-r")
# library(posteriordb)
# my_pdb <- pdb_github()
# pos <- posterior_names(my_pdb)
# print(pos)
#
# po <- "irt_2pl"
# dat <- pdb_data(po)
# code <- stan_code(po)
# #ref <- reference_posterior_draws(po)
# #df <- lapply(ref[1:4], \(x) as.data.frame(x)) |> bind_rows()
# saveRDS(dat, file='dat.RDS')
# #saveRDS(df, file='posterior.RDS')
# writeLines(code, con ='irt_2pl.stan')


setwd(here::here('models/irt_2pl'))
dat <- readRDS('dat.RDS')

library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=4)
fit.stan <- stan(file='irt_2pl.stan', data = dat, chains=4, iter=2000)
post <- as.data.frame(fit.stan)

library(shinystan)
launch_shinystan(fit.stan)
#corrplot::corrplot(cor(post), type='lower')

library(RTMB)
pars <- list(logsigma_theta=1, theta=rep(1, dat$J),
             logsigma_a=-0.8, loga=rep(.5, dat$I),
             mu_b=1, logsigma_b=1, b=rep(.5, dat$I))
func <- function(pars){
  getAll(dat, pars)
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
obj <- MakeADFun(func, pars, random=c('theta', 'loga','b'))
obj$report()$lp+func(pars)
obj$fn()

# check TMB matches Stan
p1 <- post[1,1:144]
p2 <- post[200, 1:144]
p1$sigma_a <- log(p1$sigma_a)
p1$sigma_b <- log(p1$sigma_b)
p1$sigma_theta <- log(p1$sigma_theta)
p1[103:122] <- log(p1[103:122])
p2$sigma_a <- log(p2$sigma_a)
p2$sigma_b <- log(p2$sigma_b)
p2$sigma_theta <- log(p2$sigma_theta)
p2[103:122] <- log(p2[103:122])
p1 <- as.numeric(p1)
p2 <- as.numeric(p2)
diff.tmb <- obj$report(p1)$lp -obj$report(p2)$lp
diff.stan <- post$lp__[1]-post$lp__[200]
diff.stan - diff.tmb

saveRDS(inputs, file='obj.irt_2pl.RDS')


# run longer chain and compare correlations and marginal sds
library(SparseNUTS)
obj <- MakeADFun(func, pars, random=c('theta', 'loga','b'))
fit <- sample_snuts(obj, iter=3000, warmup=500, metric='auto',
                         globals = list(dat=dat), control=list(adapt_delta=.95))
pairs(fit, pars=1:6, order='slow')
pairs(fit, pars=1:6, order='mismatch')
pairs(fit, pars=c('logsigma_a', 'loga[1]', 'loga[2]', 'loga[3]'))
pairs(fit, pars=c('logsigma_theta', 'theta[1]', 'theta[2]', 'theta[3]'))
pairs(fit, pars=c('logsigma_b', 'b[1]', 'b[2]', 'b[3]'))
post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))
hist(abs(cors[lower.tri(cors)]))

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



# how about a non-centered version?
func_nc <- function(pars){
  getAll(dat, pars)
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
obj_nc <- MakeADFun(func_nc, pars, random=c('theta', 'loga','b'))

fit_nc <- sample_snuts(obj_nc, iter=3000, warmup=500, metric='auto',
                         globals = list(dat=dat), control=list(adapt_delta=.95))
pairs(fit_nc, pars=1:6, order='slow')
pairs(fit_nc, pars=1:6, order='mismatch')
pairs(fit_nc, pars=c('logsigma_a', 'loga[1]', 'loga[2]', 'loga[3]'))
pairs(fit_nc, pars=c('logsigma_theta', 'theta[1]', 'theta[2]', 'theta[3]'))
pairs(fit_nc, pars=c('logsigma_b', 'b[1]', 'b[2]', 'b[3]'))
post <- as.data.frame(fit_nc)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))
hist(abs(cors[lower.tri(cors)]))

obj_nc$par |> length()
obj_nc$env$par |> length()
opt <- with(obj_nc, nlminb(par,fn,gr))
Q <- sdreport(obj_nc, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

minsd <- apply(post, 2, sd) |> min()
maxsd <- apply(post, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd

setwd(here::here())
