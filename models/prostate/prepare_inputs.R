
## modified from https://raw.githubusercontent.com/charlesm93/laplace_manuscript/refs/heads/master/bernoulli_logit.r

# notably converted 'prostrate_full_data.R' to dat_full.RDS

library(adnuts)
library(rstan)
rstan_options(auto_write = TRUE)
setwd(here::here('models/prostate'))



dat_full <- readRDS('dat_full.RDS')

# Use this sequence to examine a subset of the data
dat <- dat_full
dat$x <- dat$x[, (dat$d - 99):dat$d]
# If we only retain a subset of x, we need to adjust tau and d.
dat$p <- ncol(dat$x)
p0 <- 5
sigma <- with(dat, sqrt(1 / mean(y) / (1 - mean(y))))
dat$scale_global <- with(dat, p0 / (p - p0) * sigma / sqrt(n))
dat$d <- ncol(dat$x)

# verify Stan and TMB models match, take two arbitrary posterior
# vectors from Stan and ensure the lp difference is the same
# between models.
fit.stan <- rstan::stan(file='skim_logit.stan', data=dat, chains=4, cores=4,
                        iter=1000, warmup=500,
                        control=list(max_treedepth=10))
post <- as.data.frame(fit.stan)

library(TMB)
compile('prostrate.cpp')
dyn.load(dynlib('prostrate'))

xx <- as.numeric(post[1,])
nm <- names(post)
pars1 <- list()
pars1$z <- xx[grep('z', nm)]
pars1$logtau <- log(post$tau[1])
pars1$loglambda <- log(xx[grep('lambda', nm)])
pars1$logcaux <- log(post$caux[1])
pars1$logxi <- log(post$xi[1])
xx <- as.numeric(post[10,])
nm <- names(post)
pars2 <- list()
pars2$z <- xx[grep('z', nm)]
pars2$logtau <- log(post$tau[10])
pars2$loglambda <- log(xx[grep('lambda', nm)])
pars2$logcaux <- log(post$caux[10])
pars2$logxi <- log(post$xi[10])
obj1 <- MakeADFun(data=dat, parameters=pars1,
                 silent=TRUE, DLL='prostrate')
obj2 <- MakeADFun(data=dat, parameters=pars2,
                  silent=TRUE, DLL='prostrate')
post$lp__[1] - post$lp__[10]
obj1$report()$lp-obj2$report()$lp

obj <- MakeADFun(data=dat, parameters=pars1, random=c('z'),
                  silent=TRUE, DLL='prostrate')
opt <- with(obj, nlminb(par, fn, gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision


fit <- sample_sparse_tmb(obj, iter=3000, warmup=500,
                         seed=1, control=list(adapt_delta=.99))
pairs_admb(fit, pars=1:5, order='slow')
pairs_admb(fit, pars=1:5, order='mismatch')
pairs_admb(fit, pars=1:5)
plot_uncertainties(fit)
plot_Q(Q=Q)

obj$par |> length()
obj$env$par |> length() - obj$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- fit$mle$Q
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))
hist(abs(cov2cor(M)[lower.tri(M)]))

post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))


post.tmb <- as.data.frame(fit)
minsd <- apply(post.tmb, 2, sd) |> min()
maxsd <- apply(post.tmb, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd

setwd(here::here())



# # old RTMB version, fails b/c chol() errors out and crashes everything
# #library(RTMB)
# func <- function(pars){
#   RTMB::getAll(pars,dat)
#   # transformed data
#   d0 <- 5
#   slab_scale2 <- slab_scale^2
#   half_slab_df <- 0.5*slab_df
#   mu <- rep(0,n)
#   x2 <- x^2
#
#   # transformed parameters
#   tau <- exp(logtau)
#   lambda <- exp(loglambda)
#   caux <- exp(logcaux)
#   xi <- exp(logxi)
#   eta_one <- scale_global*tau
#   m_squared <- slab_scale2*caux
#   kappa_squared <- m_squared*lambda^2/(m_squared+eta_one^2*lambda^2)
#   eta_two <- eta_one^2/m_squared*xi
#   K1 <- (x%*%diag(kappa_squared)) %*% t(x)
#   K2 <- (x2%*%diag(kappa_squared)) %*% t(x2)
#   K <- 0.5*eta_two^2*(K1+1)^2 -
#     0.5*eta_two^2*K2+
#     (eta_one^2-eta_two^2)*K1 +
#     scale_icept^2 - 0.5*eta_two^2
#   diag(K) <- diag(K) + 1e-5
#   L_K <- t(chol(K))
#   func_vals <- mu+L_K%*%z
#   # jacobians
#   lp <- logtau +sum(loglambda)+logcaux + logxi +
#     # priors
#     sum(dt(lambda, nu_local, log=TRUE))+
#     dt(tau, nu_global, log=TRUE) +
#     # dinvgamma implementation directly
#     (dgamma(1/caux, shape=half_slab_df, scale=1/half_slab_df, log=TRUE)-2*log(caux))+
#     (dgamma(1/xi, shape=half_slab_df, scale=1/half_slab_df, log=TRUE)-2* log(xi))+
#     # hyperdistribution
#     sum(dnorm(z,0,1, log=TRUE)) +
#     # likelihood
#     sum(dbinom(y, size=1, prob=plogis(func_vals), log=TRUE))
#   return(-lp)
# }
#
