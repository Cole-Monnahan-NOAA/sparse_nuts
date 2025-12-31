## A quick demo and minimal reproducible example of sparse NUTS:
# -how Q is calculated by TMB using the RTMB package
# -how to use Q to transform the parameter space
# -how to pass this to Stan via StanEstimators
# -how the generalized delta method can be used for approximate
# posteriors of generated quantities

# From Monnahan, C.C., J.T. Thorson, K. Kristensen, and B.
# Carpenter (in prep). Leveraging sparsity to improve no-U-turn
# sampling efficiency for hierarchical Bayesian models

# last updated July 2025

# see https://github.com/Cole-Monnahan-NOAA/sparse_nuts/tree/main

# !! adnuts::sample_sparse_tmb has this and other functionality
# incorporated and should be used.. this is just a demo  !!

library(RTMB)
rm(list=ls())

## -------------------------------------------------------------
# Use the schools model as an example. Setup model:
dat <- list(y = c(28,  8, -3,  7, -1,  1, 18, 12),
            sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(eta=rep(1,8),mu=0, logtau=1)
f <- function(pars){
  RTMB::getAll(dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # hyperprior
    sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
    logtau                              # jacobian
  # silly generated quantity of random (eta) + fixed (tau) effects
  gq <- sum(exp(logtau)*eta[2:4])
  ADREPORT(gq)                          # request delta method SE
  REPORT(gq)                            # simple output
  return(-lp)                           # TMB uses negative lp
}
f(pars)
obj <- MakeADFun(func=f, parameters=pars, random='eta', silent=TRUE)
obj$fn()  # marginal negative log posterior
obj$gr()  # gradient of marginal -lp (fixed effects only)
# Optimize and get Q via sdreport function
opt <- with(obj, nlminb(par,fn,gr))
sdrep <- sdreport(obj, getJointPrecision=TRUE)
str(sdrep$jointPrecision)
print(round(sdrep$jointPrecision,2))
Q0 <- sdrep$jointPrecision
SparseNUTS::plot_Q(Q=Q0)

## -------------------------------------------------------------
# This code block shows how to calculate Q in R. It reproduces
# the 'sdreport' functionality in a minimal way (see equations
# 9-10) and only for RTMB models. It is meant only as a
# demonstration for the model above (schools).
library(Matrix)
q_hat <- obj$env$last.par.best # joint mode
n <- length(q_hat)
theta_hat <- opt$par          # marginal mode
r <- obj$env$random           # index of random effects
nonr <- setdiff(seq_along(q_hat), r)
# Hessian block for fixed effects using finite differences
H_Bhat <- optimHess(theta_hat, obj$fn, obj$gr) # Hessian of marginal posterior
# Hessian of random effects at joint mode using AD.
H_AA <- obj$env$spHess(q_hat, random = TRUE)
# Second derivatives of the joint posterior at the joint
# mode for the fixed:random effect elements only. Uses AD.
H_AB <- obj$env$f(q_hat, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
H_BA <- t(H_AB)
H_BB <- H_BA %*% solve(t(H_AA)) %*% H_AB + H_Bhat
Q <- Matrix(0, nrow = n, ncol=n, sparse=TRUE)
Q[r,r] <- H_AA
Q[r,nonr] <- H_AB
Q[nonr,r] <- H_BA
Q[nonr, nonr] <- H_BB
# this matches TMB's internal calculations in sdreport:
max(abs((Q-Q0))) # [1] 2.775558e-17


## -------------------------------------------------------------
# now use Q to transform the system (the sparse metric). Let 'x'
# be the original parameter space and 'y' the transformed one

# first rebuild object without the Laplace appoximation (random=NULL)
obj2 <- RTMB::MakeADFun(func=obj$env$data,
                        parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
obj2$gr() # notice this is not the marginal anymore
# Do Cholesky on permuted Q
chd <- Matrix::Cholesky(sdrep$jointPrecision, super=TRUE, perm=TRUE)
L <- as(chd, "sparseMatrix")
# Drop all numerical zeros and convert to triangular storage
L <- tril(drop0(L)) ## class(L) == "dtCMatrix"
Lt <- Matrix::t(L) ## class(Lt) == "dtCMatrix"
perm <- chd@perm + 1L            # the new order of the system
iperm <- Matrix::invPerm(perm)   # the inverse
# assume an initial value in x space, like the joint mode
x0 <- q_hat      # the joint mode
x0perm <- x0[perm]               # permuted
# initial value in the transformed, permuted space
y0perm <- as.vector(Lt %*% x0perm)
# redefine the objective and gradient functions
# equation 4 of the MS: f_q'(q')=f_q( [L'P]^{-1}q')
fn.y <- function(y)  -obj2$fn(Matrix::solve(Lt, y)[iperm])
# carefully implement this to optimize sparsity and reduce
# computations since the gradient call is the expensive part of
# NUTS
gr.y <- function(y){
  # equation 4 of the MS: g_q'(q')=g_q( [L'P]^{-1}q')[L'P]^{-1}
  x <- Matrix::solve(Lt, y)[iperm]
  -Matrix::solve(L, as.numeric(obj2$gr(x))[perm])
}
# back transform parameters
y.to.x <- function(y) as.numeric(Matrix::solve(Lt, y)[iperm])
# test them out
fn.y(y0perm)
-obj2$fn(x0)  # matches
gr.y(y0perm) # gradient at joint mode is 0 for fixed effects only

## -------------------------------------------------------------
# Now show how to link this through StanEstimators
library(StanEstimators)
fit <- stan_sample(fn=fn.y, par_inits=x0perm, grad_fun=gr.y,
                   num_chains=4, seed = 12345)
# posterior draws in transformed space
post.y <- unconstrain_draws(fit) |> as.data.frame()
# recorrelate the draws into untransformed space x
post.x <- t(apply(post.y[,2:11], 1, FUN=y.to.x))
cbind(postmean=colMeans(post.x), postmode=x0)


## -------------------------------------------------------------
# Now compare approximate (asymptotic normal) estimate of a
# generated quantity using the generalized delta method (via
# TMB::sdreport) against the posterior

# The mean and SE from the delta methd are assumed to be normal
gq.mle <- c(mean=sdrep$value, sd=sdrep$sd)
# push each posterior draw through and extract the generated
# quantity 'gq' to get posterior of gq
gq.post <- apply(post.x, 1, \(x) obj2$report(x)$gq)
hist(gq.post, freq=FALSE, breaks=30, main='Generated quantity example',
     xlab='gq')
x <- seq(min(gq.post), max(gq.post))
y <- dnorm(x, gq.mle[1], gq.mle[2])
lines(x,y, col=2, lwd=2)


# this is the package function with automates all this
mcmc <- SparseNUTS::sample_snuts(obj, globals=list(dat=dat))
plot(mcmc)
pairs(mcmc, order='slow')
