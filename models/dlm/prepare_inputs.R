
## https://james-thorson-noaa.github.io/dsem/articles/vignette.html#comparison-with-dynamic-linear-models

library(dsem)
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))

# Specify by declaring each arrow and lag
sem = "
  # Link, lag, param_name
  cprofits -> consumption, 0, a1
  cprofits -> consumption, 1, a2
  pwage -> consumption, 0, a3
  gwage -> consumption, 0, a3

  cprofits -> invest, 0, b1
  cprofits -> invest, 1, b2
  capital -> invest, 0, b3

  gnp -> pwage, 0, c2
  gnp -> pwage, 1, c3
  time -> pwage, 0, c1
"
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]
inputs <- list(sem=sem, tsdata=tsdata)

obj <- dsem(sem=inputs$sem,
           tsdata = inputs$tsdata,
           estimate_delta0 = TRUE,
            control = dsem_control(
              quiet = TRUE,
            newton_loops = 0))$obj


saveRDS(obj, file='models/dsem/obj.dlm.RDS')


fit <- sample_sparse_tmb(obj, iter=2000, warmup=1000)

pairs_admb(fit, pars=1:7, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)

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
