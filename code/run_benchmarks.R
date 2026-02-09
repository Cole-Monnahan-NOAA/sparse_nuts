source("code/startup.R")



#  benchmark of metric gradient timings on a simplified
# model, isolating rotation costs vs model + rotation gradient
# costs
metrics <- c('unit', 'diag', 'dense', 'sparse')
ndata <- floor(c(3,7, 9,13, 17, 20, 30, 40, 50,65, 80,100))
bench <- NULL; set.seed(121)
for(n in ndata){
  # number of replicates to do to stabilize estimates
  for(rr in 1:6){
    message('Simulating data and fitting model..')
    # tau not estimated since goes to 0 every once in a while and
    # breaks everything
    obj <- sim_spde_dat(n, sparse=TRUE, seed=n+rr)#, map=list(log_tau=factor(NA)))
    # this ensures good convergence, more reliable than nlminb
    opt <- TMBhelper::fit_tmb(obj, getsd=FALSE, newtonsteps = 2, loopnum=2,
                              control=list(trace=0))
    stopifnot(max(abs(obj$gr(opt$par)))<.0001)
    Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
    Qinv <- solve(Q)
    nrepars <- length(obj$env$parList()$x)
    mle <- obj$env$last.par.best
    times <- 1000#ifelse(nrepars<500, 500, 100)
    # message("Rebuilding RTMB obj without random effects...")
    for(metric in metrics){
      for(type in c('original', 'simple')){
        cat(n, rr, metric, type, '\n')
        if(type=='original'){
          # rebuilding obj each time so don't need to worry about
          # reusing x in gr(x) and messing up the timings
          obj2 <- RTMB::MakeADFun(func = obj$env$data,
                                  parameters = obj$env$parList(),
                                  map = obj$env$map, random = NULL, silent = TRUE,
                                  DLL = obj$env$DLL)
        } else {
          obj2 <- RTMB::MakeADFun(func = function(pars) sum(pars$x),
                                  parameters = obj$env$parList(),
                                  map = obj$env$map, random = NULL, silent = TRUE,
                                  DLL = obj$env$DLL)
        }
        set.seed(n)
        out <- SparseNUTS:::.rotate_posterior(metric, obj2$fn, obj2$gr, Q, Qinv, list(mle))
        size <- as.numeric(object.size(out))/1000 # Kb of memory
        # bench::mark seemed to work better but produced flat for sparse simple combo which seemed wrong.. sparse gradient calls should scale worse with dimensionality. not sure why that happened.
      #  time <- as.numeric(bench::mark(out$gr2(out$x.cur+rnorm(length(out$x.cur))/1000),
       #                   min_time=Inf, max_iterations=times, time_unit='ms')$median)
        npars <- length(mle)
        time <- as.numeric(summary(microbenchmark(grad=out$gr2(out$x.cur[[1]]),
                                       times=times, unit='ms'))$median)
          bench <- rbind(bench, data.frame(type=type,
                                         replicate=rr,
                                         metric=metric,
                                         ndata=n, nrepars=npars,
                                         size=size,
                                         time=time))
      }
    }
  }
  if(n>3){
    g <- plot.bench(bench)
    print(g)
    ggsave('plots/spde_gradient_bechmark.png', g, width=8, height=2.75, units='in')
    saveRDS(bench, 'results/bench_spde_gr.RDS')
  }
}



## Benchmark the case studies for different metrics. Used to be embedded in
## sample_snuts but split out into separate function in the package. This is
## very sensitive to things but seems to work by adding a very, very tiny random
## vector to gr() and running a ton of replicates
source('code/startup.R')
source('code/load_rtmb_objects.R')
b_m <- function(obj, model=NULL, metrics=c('unit', 'auto'), times=5000)
  SparseNUTS:::benchmark_metrics(obj=obj, model_name=model,
                                 metrics=metrics,
                                 times=times)
bench.rtmb <- bind_rows(
  b_m(obj.causal, model = 'causal'),
  b_m(obj.diamonds, model = 'diamonds'),
  b_m(obj.irt_2pl, model = 'irt_2pl'),
  b_m(obj.irt_2pl_nc, model = 'irt_2pl_nc'),
  b_m(obj.kilpisjarvi, model = 'kilpisjarvi', metrics=c('unit', 'dense')),
  b_m(obj.lynx, model = 'lynx'),
  b_m(obj.radon, model = 'radon'),
  b_m(obj.schools, model = 'schools')
)
source('code/load_tmb_objects.R')
bench.tmb <- bind_rows(
  b_m(obj.dlm, model = 'dlm'),
  b_m(obj.gp_pois_regr),
  b_m(obj.petrel, times=1000, model='petrel'),
  b_m(obj.pollock),
  b_m(obj.salamanders, model='salamanders'),
  b_m(obj.sam, model='sam'),
  b_m(obj.sdmTMB),
  b_m(obj.swallows),
  b_m(obj.wham),
  b_m(obj.wildf))

bench <- rbind(bench.rtmb, bench.tmb) |>
  group_by(model) |>
#  select(model, metric, time, npar, pct.sparsity) |>
  mutate(rel_gr=round(time/time[metric=='unit'],2)) |>
  filter(metric!='unit') |> arrange(model)

# this gets merged with other stats in process_results
write.csv(bench, file='results/bench_casestudies.csv', row.names=FALSE)

