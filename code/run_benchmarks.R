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
  mutate(rel_gr=gr/gr[metric=='unit'])
bench_wide <- select(bench, model, pct.sparsity, npar, metric, rel_gr) |>
  tidyr::pivot_wider(names_from='metric', values_from='rel_gr')
write.csv(bench_wide, file='results/bench_casestudies_table.csv', row.names=FALSE)



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
bench_spde <- readRDS('results/bench_spde.RDS') |>
  group_by(nrepars) |>
  mutate(reltime=time[model=='dense']/time)
g <- ggplot(filter(bench_spde, model=='sparse'),
            aes(nrepars, reltime, color=NULL))  + geom_line() +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(x='Number of spatial random effects',
       y='Sparse:dense\ngradient evaluation',
       color=NULL) +
  theme(legend.position.inside=c(.2,.8), legend.position='inside')
ggsave('plots/bench_spde.png', g, width=3.5, height=2.5)


# quick benchmark of metric gradient timings on a simplified model, isolating rotation costs vs model + rotation gradient costs
metrics <- c('unit', 'diag', 'dense', 'sparse-J', 'sparse', 'RTMBtape')
ndata <- floor(c(2:5,7,9, seq(10,100, len=10)))
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
  for(type in c('original', 'simple')){
    if(type=='simple'){
      obj2 <- RTMB::MakeADFun(func=function(pars) sum(pars$x), parameters=list(x=mle), silent=TRUE)
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

bench_spde_gr <- bench |> group_by(type, nrepars) |>
  mutate(reltime=time/time[metric=='unit']) |>
  mutate(type=factor(type,
                     levels=c('simple', 'original'),
                     labels = c('Rotation only', 'Rotation and model gradient')))
g <- ggplot(bench_spde_gr, aes(x=nrepars, y=reltime, color=metric)) +
  geom_line(linewidth=1, alpha=.8) + #geom_point() +
  labs(y='Time for gradient relative to unit', x='Number of random effects',
       title='Benchmark for SPDE model') +
  facet_wrap('type', nrow=2, scales='free_y') + scale_y_log10() +
  scale_x_log10()
ggsave('plots/spde_gradient_bechmark_all.png', g, width=5.5, height=5, units='in')
g <- ggplot(filter(bench_spde_gr, metric %in% c('diag', 'dense', 'sparse')),
            aes(x=nrepars, y=reltime, color=metricf(metric))) +
  geom_line(linewidth=1, alpha=.8) + #geom_point() +
  labs(y='Time for gradient relative to none', x='Number of random effects',
       title='Benchmark for SPDE model') +
  facet_wrap('type', nrow=2, scales='free_y') + scale_y_log10() +
  scale_x_log10() + labs(color='metric')
ggsave('plots/spde_gradient_bechmark.png', g, width=4.5, height=3, units='in')
saveRDS(bench_spde_gr, 'results/bench_spde_gr.RDS')
