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

map <- list(lnsigma_j=factor(c(NA,NA)))

# pad NAs to build a long time-series of simulated data below
na_mat <- matrix(NA, nrow=2048, ncol=2)
#na_mat <- matrix(NA, nrow=20, ncol=2)
data2 <- ts(rbind(data, na_mat))
fit0 = dsem( sem = sem,
             family=family,
             tsdata = data2,
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
             tsdata = data2,
             estimate_delta0 = FALSE,
             control = dsem_control(
               parameters=parameters,
               quiet = FALSE,
               run_model=TRUE,
               map=map,
               getsd = FALSE) )
parameters = fit0$obj$env$parList()
parameters$delta0_j = rep(0, ncol(data2) )

# Refit with delta0
fit = dsem( sem = sem,
            family=family,
               tsdata = data2,
               estimate_delta0 = TRUE,
               control = dsem_control( quiet=TRUE,
                                       map=map,
                                       parameters = parameters ) )


# this is the full simulated data set, which can be
datsim <- simulate(fit, nsim=1, seed=123,
                   variance='random',fill_missing = TRUE,
                   resimulate_gmrf = TRUE)[[1]]
matplot(datsim)

#test <- sample_sparse_tmb(obj, iter=1000)
# test
# pairs_admb(test, order='slow', pars=1:5)
#plot_Q(test)

options(future.globals.maxSize= 2000*1024^2) # 1000 MB limit
metrics <- c('unit', 'diag', 'dense', 'sparse')
stats <- list()
for(nyrs in c(100,250,375, 500, 750, 1500)){
  message("Starting analysis for nyrs=",nyrs)
  dattmp <- ts(datsim[1:nyrs,])
  partmp <- parameters
  partmp$x_tj <- partmp$x_tj[1:nyrs,]
  fitsim = dsem( sem = sem,
                 family=family,
                 tsdata = dattmp,
                 estimate_delta0 = TRUE,
                 control = dsem_control( quiet=TRUE,
                                         map=map,
                                         parameters = partmp ) )
  obj <- fitsim$obj
  fits <- fit_models(obj, chains=4, cores=4,
                     metrics=metrics,
                     iter=2000,
                     init='last.par.best',
                     replicates=reps, cpus=cpus,
                     model='VAR', plot=FALSE)
  nrepars <- 2*nyrs
  stats <- rbind(stats, cbind(nyrs=nyrs, nrepars=nrepars, get_stats(fits)))
  saveRDS(stats, file='results/VAR_stats.RDS')
  g <- ggplot(stats, aes(nrepars, eff, group=interaction(metric,replicate), color=metric)) +
    geom_point(alpha=.5) +
    geom_line() +
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
