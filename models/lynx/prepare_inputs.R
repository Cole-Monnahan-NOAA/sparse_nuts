# Lynx hare dynamics modeled as ODE with RK4 solver. Started by
# J. Thorson, modified from there.
library(RTMB)

# local copy of ecostate::rk4sys
rk4sys <- function (f, a, b, y0, n, Pars, ...){
  m <- length(y0)
  h <- (b - a)/n
  x <- seq(a + h, b, by = h)
  y <- matrix(0, nrow = n, ncol = m)
  k1 <- h * f(a, y0, Pars, ...)
  k2 <- h * f(a + h/2, y0 + k1/2, Pars, ...)
  k3 <- h * f(a + h/2, y0 + k2/2, Pars, ...)
  k4 <- h * f(a + h, y0 + k3, Pars, ...)
  y[1, ] <- y0 + k1/6 + k2/3 + k3/3 + k4/6
  for (i in seq_len(n - 1)) {
    k1 <- h * f(x[i], y[i, ], Pars, ...)
    k2 <- h * f(x[i] + h/2, y[i, ] + k1/2, Pars, ...)
    k3 <- h * f(x[i] + h/2, y[i, ] + k2/2, Pars, ...)
    k4 <- h * f(x[i] + h, y[i, ] + k3, Pars, ...)
    y[i + 1, ] <- y[i, ] + k1/6 + k2/3 + k3/3 + k4/6
  }
  x <- c(a, x)
  y <- rbind(y0, y)
  return(list(x = x, y = y))
}
f <- identity
dzdt <- function(Time, State, Pars){
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  d1_dt = f(Pars$ln_alpha) * State[1] - f(Pars$ln_beta) * State[1] * State[2]
  d2_dt = -1 * f(Pars$ln_gamma) * State[2] + f(Pars$ln_delta) * State[1] * State[2]
  return( c(d1_dt, d2_dt) )
}

get_nll <- function(pars){
  getAll(pars)
  zhat_ti <- array(0, dim = dim(dat) )
  # Necessary in packages
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  # Initialize
  zhat_ti[1,] <- ln_z0
  for( t_index in 2:nrow(dat) ){
    proj <- rk4sys(f=dzdt, a=0, b=1, n=10,
                   Pars=pars, y0=z_ti[t_index-1,])
    zhat_ti[t_index,] <- proj$y[nrow(proj$y),]
  }
  lp <- ln_sigma + ln_tau + # jacobians
    sum(dnorm(dat, mean=z_ti, sd=exp(ln_sigma), log=TRUE)) +
    sum(dnorm(z_ti, mean=zhat_ti, sd=exp(ln_tau), log=TRUE)) +
    # informative prior
    dgamma(x=exp(ln_sigma), scale=.03, shape=9.5, log=TRUE)
  ADREPORT( zhat_ti )
  REPORT(lp)
  return(-lp)
}

pars <- list(ln_alpha=1, ln_beta=1, ln_gamma=1,
             ln_delta=1, ln_z0=c(0.1,0.1),
             ln_sigma=1, ln_tau=1,
             z_ti = dat) #ifelse( is.na(dat), rep(1,nrow(dat))%o%colMeans(dat,na.rm=TRUE), dat
# https://jmahaffy.sdsu.edu/courses/f09/math636/lectures/lotka/qualde2.html
hares <- c(30, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22, 25.4,
           27.1, 40.3, 57, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7)
lynx <- c(4, 6.1, 9.8, 35.2, 59.4, 41.7, 19, 13, 8.3, 9.1, 7.4,
          8, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6)
dat <- cbind(X=hares/10, Y=lynx/10)

get_nll(pars)
obj <- MakeADFun( func = get_nll, silent=TRUE,
                 parameters = pars,
                 random = "z_ti" )
opt <- nlminb( obj$par, obj$fn, obj$gr, # hessian = obj$he,
              control = list(iter.max = 1e4, eval.max = 1e4))



library(adnuts)
globals <- list(dat=dat, rk4sys=rk4sys, dzdt=dzdt, f=f)
fit <- sample_sparse_tmb(obj, seed=1, globals=globals,
                         iter=3000, warmup=500,
                         control=list(adapt_delta=.99))

pairs_admb(fit, order='mismatch', pars=1:5)
pairs_admb(fit, order='slow', pars=1:5)
pairs_admb(fit, pars=c('ln_sigma', 'ln_tau', 'z_ti[1]', 'z_ti[2]','z_ti[3]'))

post <- as.data.frame(fit)
hist(exp(post$ln_sigma), freq=FALSE)
x <- seq(.0001,1, len=1000)
y <- dgamma(x=x,shape=9.5, scale=.03)
lines(x,y)

obj$par |> length()
obj$env$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

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




# # old explorations of likelihod profiles
# ## 2d profile over the sigmas
# prof <- TMB::tmbprofile(obj, name='ln_sigma')
# plot(prof)
#
# df <- expand.grid(ln_sigma=seq(-5,-.5, len=30),
#                   ln_tau=seq(-2,-.5, len=30),
#                   include_prior=c(TRUE,FALSE),
#                   nll=NA)
# for(i in 1:nrow(df)){
#   parlist$ln_sigma <- df$ln_sigma[i]
#   parlist$ln_tau <- df$ln_tau[i]
#   include_prior <- df$include_prior[i]
#   objtmp <- MakeADFun( func = get_nll, silent=TRUE,
#                        parameters = parlist,
#                        map=list(ln_tau=factor(NA), ln_sigma=factor(NA)),
#                        random = "z_ti" )
#   opttmp <- with(objtmp, nlminb(par, fn, gr))
#   df$nll[i] <- opttmp$objective
# }
# library(dplyr)
# df <- group_by(df, include_prior) |>
#   mutate(delta_nll=nll-min(nll))
# #df$delta_nll <- pmin(df$delta_nll,4)
#
# library(ggplot2)
# ggplot(subset(df, delta_nll<4), aes(ln_sigma, ln_tau, color=delta_nll))  +
#   geom_point() + facet_wrap('include_prior')
