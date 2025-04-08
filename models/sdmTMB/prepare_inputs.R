setwd(here::here('models/sdmTMB/'))


library(sdmTMB)
## Example model
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
fit0 <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

## refit with ln_kappa mapped off since inestimable as is
obj <- fit0$tmb_obj
pars2 <- fit0$parlist
map2 <- fit0$tmb_map
map2$ln_kappa <- factor(c(NA,NA))
obj <- TMB::MakeADFun(data=obj$env$data, parameters=pars2,
                             map=map2, random=fit0$tmb_random, silent=TRUE, DLL=obj$env$DLL)


saveRDS(obj, file='obj.sdmTMB.RDS')


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
