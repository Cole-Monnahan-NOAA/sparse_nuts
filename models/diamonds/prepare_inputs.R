# ## remotes::install_github("stan-dev/posteriordb-r")
# library(posteriordb)
# my_pdb <- pdb_github()
# po <- 'diamonds'
# dat <- pdb_data(po)
# code <- stan_code(po)
# writeLines(code, con='models/diamonds/diamonds.stan')
# saveRDS(dat, file='models/diamonds/dat.RDS')

library(RTMB)

setwd(here::here('models/diamonds'))
diamonds_dat <- readRDS('dat.RDS')

library(rstan)
fit.stan <- stan('diamonds.stan', model_name='diamonds',
                 data=diamonds_dat, cores=4)


# Transform data before passing to TMB
Kc <- diamonds_dat$K - 1;
Xc <- matrix(NA, nrow=diamonds_dat$N, ncol=Kc) # centered version of X without an intercept
means_X <- rep(NA, Kc) # column means of X before centering
for(i in 2 : diamonds_dat$K) {
  means_X[i - 1] = mean(diamonds_dat$X[,i]);
  Xc[,i - 1] = diamonds_dat$X[,i] - means_X[i - 1];
}
diamonds_dat$means_X <- means_X
diamonds_dat$Xc <- Xc
pars <- list(b=rep(0, Kc),## population-level effects
             Intercept=1,   ## temporary intercept for centered predictors
             logsigma=0) ## residual SD

## define the (negative) posterior density as a function in R
f <- function(pars){
  getAll(pars,diamonds_dat)
  ## transformed parameters
  sigma <- exp(logsigma)
  lp <- logsigma + ## Jacobian
    ## priors
    sum(dnorm(b, 0, 1, log=TRUE)) +
    sum(dt((Intercept-8)/10, 3,log=TRUE) )+ # dropped constants
    sum(dt(sigma/10,df=3, log=TRUE))
  ## likelihood including all constants
  mu_hat <- Intercept + as.numeric(Xc%*%b)
  if (!prior_only)  lp <- lp+sum(dnorm(Y, mean=mu_hat, sd=sigma, log=TRUE))
  ## generated quantities
  REPORT(mu_hat)
  b_Intercept <- Intercept - sum(means_X*b)
  REPORT(b_Intercept)
  REPORT(lp)
  return(-lp) # TMB expects negative log posterior
}
obj <- RTMB::MakeADFun(f, pars, random='b', silent=TRUE)



# check TMB matches Stan
p1 <- post[1,-390]
p2 <- post[1000,-390]
diff.stan <- post$lp__[1]-post$lp[1000]
p1$sigma <- log(p1$sigma)
p2$sigma <- log(p2$sigma)
diff.tmb <- obj$report(as.numeric(p1))$lp-
  obj$report(as.numeric(p2))$lp
diff.stan - diff.tmb

## run longer chains
library(adnuts)
library(StanEstimators)
fit <- sample_sparse_tmb(obj, iter=2000, warmup=200, chains=4,
                         cores=4,
                         globals=list(diamonds_dat=diamonds_dat))

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



pairs_admb(fit, pars=1:6, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)

setwd(here::here())
