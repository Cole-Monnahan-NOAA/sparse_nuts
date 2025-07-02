# An RTMB port of the "irt_2pl" model in posteriordb-r data base

# library(posteriordb)
# my_pdb <- pdb_github()
# code <- stan_code('irt_2pl')
# print(code)

dat <- readRDS('dat.RDS')
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

opt <- with(obj, nlminb(par,fn,gr))
print(sdreport(obj))
