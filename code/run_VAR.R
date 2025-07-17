library(dsem)
data(isle_royale)
data = ts( log(isle_royale[,2:3]), start=1959)
family=c('normal', 'normal')

sem = "
  # Link, lag, param_name
  wolves -> wolves, 1, arW
  moose -> wolves, 1, MtoW
  wolves -> moose, 1, WtoM
  moose -> moose, 1, arM
"

# initial fit to build model with ln_sigma fixed
map <- list(lnsigma_j=factor(c(NA,NA)))
fit0 = dsem( sem = sem,
             family=family,
             tsdata = data,
             estimate_delta0 = FALSE,
             control = dsem_control(
               quiet = FALSE,
               run_model=FALSE,
               getsd = FALSE) )
parameters = fit0$obj$env$parList()
parameters$lnsigma_j <- rep(log(.1),2)

# initial fit to get good starting values
fit0 = dsem( sem = sem,
             family=family,
             tsdata = data,
             estimate_delta0 = FALSE,
             control = dsem_control(
               parameters=parameters,
               quiet = FALSE,
               run_model=TRUE,
               map=map,
               getsd = FALSE) )
parameters = fit0$obj$env$parList()
parameters$delta0_j = rep(0, ncol(data) )

# Refit with delta0
fit = dsem( sem = sem,
            family=family,
            tsdata = data,
            estimate_delta0 = TRUE,
            control = dsem_control( quiet=TRUE,
                                    map=map,
                                    parameters = parameters ) )
parameters = fit$obj$env$parList()

options(future.globals.maxSize= 5000*1024^2) # 1000 MB limit
metrics <- c('unit', 'auto')
stats <- list()
for(nyrs in  c(100,250,375, 500, 750, 1500)){
  message("Starting analysis for nyrs=",nyrs)
  # pad NAs to build a long time-series of simulated data below
  na_mat <- matrix(NA, nrow=nyrs-61, ncol=2)
  dattmp <- ts(rbind(data, na_mat))
  partmp <- parameters
  partmp$x_tj <- matrix(0, ncol=2, nrow=nrow(dattmp))
  message("Simulating data..")
  # get dimensions right
  fitsim = dsem( sem = sem,
                 family=family,
                 tsdata = dattmp,
                 estimate_delta0 = TRUE,
                 control = dsem_control( quiet=TRUE,
                                         map=map,
                                         parameters = partmp ) )
  # simulate a whole new data set
  datsim <- simulate(fitsim, nsim=1, seed=nyrs,
                     variance='random',fill_missing = TRUE,
                     resimulate_gmrf = TRUE)[[1]]
  # rebuild to get the data into the obj for MCMC sampling
  fitsim <- dsem(sem = sem, family=family, tsdata = datsim,
                 estimate_delta0 = TRUE,
                 control = dsem_control(quiet=TRUE,
                                         map=map,
                                         parameters = partmp ) )
  obj <- fitsim$obj
  fits <- fit_models(obj, chains=4, cores=4,
                     metrics=metrics,
                     iter=2000,
                     init='last.par.best',
                     control=list(adapt_delta=.95,max_treedepth=1),
                     replicates=reps, cpus=cpus,
                     model='VAR', plot=FALSE)
  nrepars <- 2*nyrs
  stats <- rbind(stats, cbind(nyrs=nyrs, nrepars=nrepars, get_stats(fits)))
  saveRDS(stats, file='results/VAR_stats.RDS')
  g <- ggplot(stats, aes(nrepars, eff, group=interaction(metric,replicate), color=metric)) +
    geom_point(alpha=.5) +
    geom_line() + scale_x_log10()+
    scale_y_log10() + labs(x='# of random effects', y='Efficiency (ESS/t)')
  ggsave('plots/VAR_stats.png', g, width=7, height=5)
  print(g)
}

# M <- fits[[1]]$mle$Qinv
# corrplot::corrplot(cov2cor(M), type='upper')
# max(abs(cov2cor(M)[lower.tri(M)]))
# hist(abs(cov2cor(M)[lower.tri(M)]))
# minsd <- min(sqrt(diag(M)))
# maxsd <- max(sqrt(diag(M)))
# maxsd/minsd
