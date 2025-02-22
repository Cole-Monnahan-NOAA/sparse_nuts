source("code/startup.R")
require(microbenchmark)
library(Matrix)
#detach(package:TMB)
library(RTMB)

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

ndata <- seq(3,100, by=5)
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
    if(type=='simple') obj2$gr <- function(x) t(x)

    runit <- adnuts:::.rotate_posterior('unit', obj2$fn, obj2$gr, Q, Qinv, mle)
    rdiag <- adnuts:::.rotate_posterior('diag', obj2$fn, obj2$gr, Q, Qinv, mle)
    rdense <- adnuts:::.rotate_posterior('dense', obj2$fn, obj2$gr, Q, Qinv, mle)
    rsparse <- adnuts:::.rotate_posterior('sparse', obj2$fn, obj2$gr, Q, Qinv, mle)
    # rauto <- adnuts:::.rotate_posterior('auto', obj2$fn, obj2$gr, Q, Qinv, mle)


    xx <- microbenchmark(runit$gr2(runit$x.cur),
                         rdiag$gr2(rdiag$x.cur),
                         rdense$gr2(rdense$x.cur),
                         rsparse$gr2(rsparse$x.cur),
                         times=200, unit='s')
    bench <- rbind(bench,
                   data.frame(type=type,
                              metric=c('unit', 'diag', 'dense', 'sparse'),
                              ndata=n, nrepars=nrepars,
                              time=summary(xx)$median))

  }
  bench2 <- bench |> group_by(type, nrepars) |>
    mutate(reltime=time/time[metric=='unit'])

  g <- ggplot(bench2, aes(x=nrepars, y=time, color=metric)) +
    geom_line() + geom_point() +
    facet_wrap('type')+
    scale_y_log10() + scale_x_log10()
  print(g)
}


bench2 <- bench |> group_by(type, nrepars) |>
  mutate(reltime=time/time[metric=='unit'])
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
ggsave('plots/spde_gradient_bechmark.png', g, width=4, height=5, units='in')
saveRDS(bench, 'results/bench_spde_rotation.RDS')

# FYI I've modified my code to default to an "auto" metric which uses an algorithm to determine which one to use before starting sampling, including a test for which gradient calc is more expensive. As such it tends to select "dense" at lower dimensions and "sparse" at higher when there are meaningful correlations (and select 'diag' when not). 
# Here's the current benchmark output for a SPDE Poisson simulation model:It's nice to see the blue line flatten out which is the main benefit of using a sparse metric. We can decorrelate the posterior with very little overhead. If we could get the earlier part to befaster with this way that would be great.
