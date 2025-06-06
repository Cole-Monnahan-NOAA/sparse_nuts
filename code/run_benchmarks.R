source("code/startup.R")
require(microbenchmark)
library(Matrix)
detach(package:TMB)
library(RTMB)

## benchmark the case studies for different metrics
source('code/load_rtmb_objects.R')
b_m <- benchmark_metrics
bench.rtmb <- bind_rows(
  b_m(obj.causal, model = 'causal'),
  b_m(obj.diamonds, model = 'diamonds'),
  b_m(obj.irt_2pl, model = 'irt_2pl'),
  #b_m(obj.kilpisjarvi, model = 'kilpisjarvi')
  b_m(obj.lynx, model = 'lynx'),
  b_m(obj.radon, model = 'radon'),
  b_m(obj.schools, model = 'schools')
)
source('code/load_tmb_objects.R')
objs <- list(obj.gp_pois_regr, obj.petrel, obj.pollock, obj.salamanders, obj.sam, obj.sdmTMB, obj.simple, obj.swallows, obj.wham, obj.wildf)
bench.tmb <-  objs |> lapply(b_m) |> bind_rows()
bench <- rbind(bench.rtmb, bench.tmb)
saveRDS(bench, file='results/bench_casestudies.RDS')

bench <- readRDS(file='results/bench_casestudies.RDS') |>
  group_by(model) |>
  mutate(rel_gr=gr/gr[metric=='sparse'])
filter(bench, metric %in% c('RTMBtape', 'sparse', 'sparse-J'))
ggplot(bench, aes(x=model, y=rel_gr, group=metric, color=metric)) + geom_line() + scale_y_log10()
ggplot(bench, aes(x=npar, y=gr, shape=metric, color=model)) + geom_point() +
  scale_y_log10() + scale_x_log10()
ggplot(bench, aes(x=metric, y=rel_gr, group=model,color=model)) + geom_line() +
  scale_y_log10()

## benchmark gradient calcs externally to show the impact of
## sparsity on the joint grad
bench.sparse <- bench.dense <- list()
ndata <- floor(sqrt(10^(1:3)))
ndata <- floor(sqrt(2^(4:11)))
for(n in ndata){
  print(n)
  obj <- sim_spde_dat(n, sparse=TRUE)
  nrepars <- length(obj$env$parList()$x)
  obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                          map=obj$env$map,
                          random=NULL, silent=TRUE,
                          DLL=obj$env$DLL)
  bench <- microbenchmark(obj2$gr(), times=50, unit='s')
  bench.sparse <- rbind(bench.sparse,
                        data.frame(model='sparse',
                                   ndata=n, nrepars=nrepars,
                                   time=summary(bench)$median[1]))
  obj <- sim_spde_dat(n, sparse=FALSE)
  obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                          map=obj$env$map,
                          random=NULL, silent=TRUE,
                          DLL=obj$env$DLL)
  bench <- microbenchmark(obj2$gr(), times=50, unit='s')
  bench.dense <- rbind(bench.dense,
                       data.frame(model='dense',
                                  ndata=n, nrepars=nrepars,
                                  time=summary(bench)$median[1]))
}
bench <- rbind(bench.sparse, bench.dense)
saveRDS(bench, 'results/bench_spde.RDS')


# quick benchmark of metric gradient timings
metrics <- c('RTMBtape-perm', 'unit', 'diag', 'sparse-Jnoperm',  'dense', 'sparse-perm-inefficient',
             'sparse-perm-efficient')
metrics <- c('unit', 'diag', 'dense', 'sparse-Jnoperm', 'sparse-perm-efficient', 'RTMBtape-perm')
ndata <- seq(2,75, by=8)
bench <- NULL
for(n in ndata){
  print(n)
  obj <- sim_spde_dat(n, sparse=TRUE)
  opt <- with(obj, nlminb(par, fn, gr))
  Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
  Qinv <- solve(Q)
  nrepars <- length(obj$env$parList()$x)
  mle <- obj$env$last.par.best
  # message("Rebuilding RTMB obj without random effects...")
  obj2 <- RTMB::MakeADFun(func = obj$env$data, parameters = obj$env$parList(),
                          map = obj$env$map, random = NULL, silent = TRUE,
                          DLL = obj$env$DLL)
  for(type in c('original', 'simple')[2]){
    if(type=='simple'){
      obj2 <- RTMB::MakeADFun(func=function(pars) sum(pars$x^2), parameters=list(x=mle), silent=TRUE)
    }
    for(metric in metrics){
      out <- adnuts:::.rotate_posterior(metric, obj2$fn, obj2$gr, Q, Qinv, mle, obj=obj2)
      xx <- microbenchmark(grad=out$gr2(out$x.cur+rnorm(length(mle),sd=.1)),
                           times=200, unit='ms')
      bench <- rbind(bench, data.frame(type=type,
                                       metric=metric,
                                       ndata=n, nrepars=nrepars, time=summary(xx)$median))
    }

  }
  if(n>2){
    bench2 <- bench |> group_by(type, nrepars) |> filter(type=='simple') |>
      mutate(reltime=time/time[metric=='unit'])
    g <- ggplot(bench2, aes(x=nrepars, y=time, color=metric)) +
      geom_line() + geom_point() +
      facet_wrap('type')+
      scale_y_log10() + scale_x_log10()
    print(g)
  }
}

bench2 <- mutate(bench2, metric=metricf(metric),
                 type=factor(type,
                             levels=c('simple', 'original'),
                             labels = c('Rotation only', 'Rotation and model gradient')))
g <- ggplot(bench2, aes(x=nrepars, y=reltime, color=metric)) +
  geom_line(linewidth=1, alpha=.8) + #geom_point() +
  labs(y='Time for gradient relative to unit', x='Number of random effects',
       title='Benchmark for SPDE model') +
  facet_wrap('type', nrow=2, scales='free_y') + scale_y_log10() +
  scale_x_log10()
g
ggsave('plots/spde_gradient_bechmark.png', g, width=4, height=5, units='in')
saveRDS(bench, 'results/bench_spde_rotation.RDS')

# FYI I've modified my code to default to an "auto" metric which uses an algorithm to determine which one to use before starting sampling, including a test for which gradient calc is more expensive. As such it tends to select "dense" at lower dimensions and "sparse" at higher when there are meaningful correlations (and select 'diag' when not). 
# Here's the current benchmark output for a SPDE Poisson simulation model:It's nice to see the blue line flatten out which is the main benefit of using a sparse metric. We can decorrelate the posterior with very little overhead. If we could get the earlier part to befaster with this way that would be great.
