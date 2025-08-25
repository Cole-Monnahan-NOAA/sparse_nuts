

library(glmmTMB)

setwd(here('models/salamanders'))

data(Salamanders)
mod = glmmTMB(count~spp * mined + (1|site), Salamanders, family="nbinom2")

## refit with ln_kappa mapped off since inestimable as is
obj <- mod$obj
saveRDS(obj, file='obj.salamanders.RDS')


fit <- sample_sparse_tmb(obj, iter=3000, warmup=500)
pairs_admb(fit, pars=1:5, order='slow')

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
