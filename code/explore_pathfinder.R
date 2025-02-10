library(StanEstimators)
library(RTMB)
library(shinystan)
library(transport)

dat <- list(J = 8,
            y = c(28,  8, -3,  7, -1,  1, 18, 12),
            sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(mu=0, logtau=1, eta=rep(1,8))
map <- NULL # estimate tau
#map <- list(logtau=factor(NA)) # do not estimate tau

f <- function(pars){
  RTMB::getAll(dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # hyperprior
    sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
    logtau                              # jacobian
  return(-lp)                           # TMB uses negative lp
}
# joint object for Stan
f(pars)
obj <- MakeADFun(func=f, parameters=pars, map=map, silent=TRUE)

fn <- function(x) -obj$fn(x)
grad_fun <- function(x) -obj$gr(x)
# marginal object for TMB
obj2 <- MakeADFun(func=f, parameters=pars, map=map,
                  random='eta', silent=TRUE)
#optimize and get Q and M=Q^(-1)
opt <- with(obj2, nlminb(par,fn,gr))
sdrep <- sdreport(obj2, getJointPrecision=TRUE)
Q <- sdrep$jointPrecision
M <- solve(Q) |> as.matrix()
ind <- seq_along(obj$par)
parnames <- names(obj$par)
parnames <- sapply(unique(names(obj$par)),
                   \(x) {
                     temp <- parnames[parnames == x]
                     if (length(temp) > 1)
                       paste0(temp, "[", 1:length(temp), "]")
                     else temp
                   }) |> unlist()

# get MCMC posterior draws for reference
fit <- stan_sample(fn=fn, grad_fun=grad_fun,
                  par_inits=obj$env$par,
                  num_chains = 1, thin=1, refresh = 1000,
                  num_samples=10000, num_warmup=1000)

# run pathfinder with default of 1000 draws.. works better if I
# cheat and start PF from the TMB mode but here starting from
# same place as TMB optimization (obj$par)
pf <- stan_pathfinder(fn=fn, grad_fun=grad_fun,
                      par_inits=obj$par)

# plot posterior (black), pathfinder (red), and Q draws (green)
library(dplyr)
draws.post <- fit@draws |> as.data.frame() |> select(1+ind) |>
  setNames(parnames) |> mutate(type='posterior')
draws.pf <- pf@draws |> as.data.frame() |> select(2+ind) |>
  mutate(type='pathfinder')
draws.Q <- mvtnorm::rmvnorm(n=1000, mean=obj2$env$last.par.best, sigma=M) |>
  as.data.frame() |> cbind(type='Qrandom')
names(draws.pf) <- names(draws.Q) <- names(draws.post)
draws2 <- rbind(draws.pf, draws.Q) |> slice_sample(prop = 1)
draws <- rbind(draws.post, draws2) |>
  mutate(col=as.numeric(factor(type, levels=c('posterior', 'pathfinder', 'Qrandom'))),
         cex=ifelse(type=='posterior', .5, 1))
pairs(draws[,1:4], pch=1, cex=draws$cex, col=draws$col, upper.panel = NULL)

## calculate wasserstein distance
library(transport)
a.pf = wpp(draws.pf[,ind], mass = rep(1 / nrow(draws.pf), nrow(draws.pf)))
a.Q = wpp(draws.Q[,ind], mass = rep(1 / nrow(draws.Q), nrow(draws.Q)))
b = wpp(draws.post[,ind], mass = rep(1 / nrow(draws.post), nrow(draws.post)))
w.pf <- wasserstein(a=a.pf, b=b, p = 1)
w.Q <- wasserstein(a=a.Q, b=b, p = 1)

w.pf
w.Q
