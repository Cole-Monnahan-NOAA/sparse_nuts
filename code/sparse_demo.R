# A quick demo of how Q is calculated by TMB using the RTMB package


library(RTMB)
rm(list=ls())
# use the schools model as an example
dat <- list(J = 8,
            y = c(28,  8, -3,  7, -1,  1, 18, 12),
            sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(eta=rep(1,8),mu=0, logtau=1)
f <- function(pars){
  RTMB::getAll(dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # hyperprior
    sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
    logtau                              # jacobian
  return(-lp)                           # TMB uses negative lp
}
obj <- MakeADFun(func=f, parameters=pars, random='eta', silent=TRUE)
#optimize and get Q and M=Q^(-1)
opt <- with(obj, nlminb(par,fn,gr))
# ask TMB to return Q from sdreport
sdrep <- sdreport(obj, getJointPrecision=TRUE)



# This code block reproduces the 'sdreport' functionality in a
# minimal way and only for RTMB models. It is meant only as a
# demonstration for the model above (schools).
library(Matrix)
par <- obj$env$last.par.best
par.fixed <- opt$par
r <- obj$env$random
nonr <- setdiff(seq_along(par), r)
hessian.fixed <- optimHess(par.fixed, obj$fn, obj$gr)
hessian.random <- obj$env$spHess(par, random = TRUE)
#obj$env$f(par, order = 0, type = "ADGrad") ## NOTE_2: ADGrad forward sweep now initialized !
tmp <- obj$env$f(par, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
A <- solve(hessian.random, tmp)
G <- as.matrix(hessian.random %*% A)
M1 <- cbind2(hessian.random, G)
M2 <- cbind2(t(G), as.matrix(t(A) %*% G) + hessian.fixed)
M <- rbind2(M1, M2)
M <- forceSymmetric(M, uplo = "L")
p <- invPerm(c(r, (1:length(par))[-r]))
# the joint precision Q
Q <- M[p, p]
class(Q) # dsCMatrix from Matrix package
Matrix::image(Q)
# this matches sdreport
max(abs((Q-sdrep$jointPrecision)))


# now use Q to transform the system (the sparse metric). Let 'x'
# be the original parameter space and 'y' the transformed one

# first rebuild object without the Laplace appoximation turned on
obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
# Do Cholesky on permuted Q
chd <- Matrix::Cholesky(sdrep$jointPrecision, super=TRUE, perm=TRUE)
L <- as(chd, "sparseMatrix")
# Drop all numerical zeros and convert to triangular storage
L <- tril(drop0(L)) ## class(L) == "dtCMatrix"
Lt <- Matrix::t(L) ## class(Lt) == "dtCMatrix"
perm <- chd@perm + 1L            # the new order of the system
iperm <- Matrix::invPerm(perm)   # the inverse
# assume an initial value in x space, like the joint mode
x0 <- obj$env$last.par.best      # the joint mode
x0perm <- x0[perm]               # permuted
# initial value in the transformed, permuted space
y0perm <- as.vector(Lt %*% x0perm)
# redefine the objective and gradient functions
fn.y <- function(y)  -obj2$fn(Matrix::solve(Lt, y)[iperm])
# carefully implement this to optimize sparsity and reduce
# computations since the gradient call is the expensive part of
# NUTS
gr.y <- function(y){
  x <- Matrix::solve(Lt, y)[iperm]
  -Matrix::solve(L, as.numeric(obj2$gr(x))[perm])
}
# back transform parameters
y.to.x <- function(y)   as.numeric(Matrix::solve(Lt, y)[iperm])

fn.y(y0perm)
obj2$fn(x0)  # matches but for sign
gr.y(y0perm) # gradient at joint mode is 0 for fixed effects only


# Now show how to link this through StanEstimators
library(StanEstimators)
fit <- stan_sample(fn=fn.y, par_inits=x0perm, grad_fun=gr.y,
                   num_chains=1)
# posterior draws in transformed space
post.y <- unconstrain_draws(fit) |> as.data.frame()
# recorrelate them
post.x <- t(apply(post.y[,2:11], 1, FUN=y.to.x))
cbind(postmean=colMeans(post.x), postmode=x0)
