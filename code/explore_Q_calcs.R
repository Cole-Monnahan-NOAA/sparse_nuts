library(StanEstimators)
library(RTMB)
library(shinystan)

dat <- list(J = 8,
            y = c(28,  8, -3,  7, -1,  1, 18, 12),
            sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(eta=rep(1,8),mu=0, logtau=1)
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
obj <- MakeADFun(func=f, parameters=pars, map=map,
                  random='eta', silent=TRUE)
#optimize and get Q and M=Q^(-1)
opt <- with(obj, nlminb(par,fn,gr))
sdrep <- sdreport(obj, getJointPrecision=TRUE)
Q <- sdrep$jointPrecision |> as.matrix()

cov.all <- matrix(0, nrow=10, ncol=10)
Huu <- numDeriv::hessian(f=obj$env$f, x=obj$env$last.par.best)[1:8, 1:8]
ff <- function(theta){
  obj$fn(theta)
  obj$env$last.par
}
J <- numDeriv::jacobian(ff, x=opt$par)
V <- sdreport(obj)$cov.fixed
cov.all[1:8, 1:8] <- solve(Huu)
cov.all <- cov.all + J %*% V %*% t(J)
cov.all - solve(Q)


## trying to steal pieces from sdreport to simplify down and get
## the sparse calcs..

# par <- obj$env$last.par.best
# par.fixed <- opt$par
# obj2 <- MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun",
#                   ADreport = TRUE, DLL = obj$env$DLL,
#                   silent =obj$env$silent)
# r <- obj$env$random
# gradient.fixed <- obj$gr(par.fixed)
# hessian.fixed <- optimHess(par.fixed, obj$fn, obj$gr)
# pdHess <- !is.character(try(chol(hessian.fixed), silent = TRUE))
# Vtheta <- try(solve(hessian.fixed), silent = TRUE)
# hessian.random <- obj$env$spHess(par, random = TRUE)
# L <- obj$env$L.created.by.newton
# phi <- try(obj2$fn(par), silent = TRUE)
# # if (is.character(phi) | length(phi) == 0) {
# #   phi <- numeric(0)
# # }
# ADGradForward0Initialized <- FALSE
# ADGradForward0Initialize <- function() {
#   obj$env$f(par, order = 0, type = "ADGrad")
#   ADGradForward0Initialized <<- TRUE
# }
# if (is.null(chunk)) {
#   Dphi <- obj2$gr(par)
# } else {
#   w <- rep(0, length(phi))
#   phiDeriv <- function(i) {
#     w[i] <- 1
#     obj2$env$f(par, order = 1, rangeweight = w,
#                doforward = 0)
#   }
#   Dphi <- t(sapply(chunk, phiDeriv))
#   phi <- phi[chunk]
# }
# if (!is.null(r)) {
#   Dphi.random <- Dphi[, r, drop = FALSE]
#   Dphi.fixed <- Dphi[, -r, drop = FALSE]
#   if (all(Dphi.random == 0)) {
#     simpleCase <- TRUE
#     Dphi <- Dphi.fixed
#   }
# }
#
# if (simpleCase) {
#   if (length(phi) > 0) {
#     cov <- Dphi %*% Vtheta %*% t(Dphi)
#   }
#   else cov <- matrix(, 0, 0)
# }
# else {
#   tmp <- solve(hessian.random, t(Dphi.random))
#   tmp <- as.matrix(tmp)
#   term1 <- Dphi.random %*% tmp
#   if (ignore.parm.uncertainty) {
#     term2 <- 0
#   }
#   else {
#     f <- obj$env$f
#     w <- rep(0, length(par))
#     if (!ADGradForward0Initialized)
#       ADGradForward0Initialize()
#     reverse.sweep <- function(i) {
#       w[r] <- tmp[, i]
#       -f(par, order = 1, type = "ADGrad", rangeweight = w,
#          doforward = 0)[-r]
#     }
#     A <- t(do.call("cbind", lapply(seq_along(phi),
#                                    reverse.sweep))) + Dphi.fixed
#     term2 <- A %*% (Vtheta %*% t(A))
#   }
#   cov <- term1 + term2
# }
# cov
#
# ans <- list(value = phi, sd = sd, cov = cov, par.fixed = par.fixed,
#             cov.fixed = Vtheta, pdHess = pdHess, gradient.fixed = gradient.fixed)
#
#
#
# hessian.random <- obj$env$spHess(par, random = TRUE)
# G <- hessian.random %*% A
# G <- as.matrix(G)
# M1 <- cbind2(hessian.random, G)
# M2 <- cbind2(t(G), as.matrix(t(A) %*% G) +
#                hessian.fixed)
# M <- rbind2(M1, M2)
# M <- forceSymmetric(M, uplo = "L")
# dn <- c(names(par)[r], names(par[-r]))
# dimnames(M) <- list(dn, dn)
# p <- invPerm(c(r, (1:length(par))[-r]))
# ans$jointPrecision <- M[p, p]
