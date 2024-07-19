


stats <- list()
for(ndata in c(10,15,20,25,30,35)){
  message("Starting analysis for ndata=",ndata)
  obj <- sim_spde_dat(ndata, TRUE, list(log_tau=factor(NA)))
  fits <- fit_models(obj, iter=1000, warmup=200, chains=4, cores=4,
                     init='random', replicates=reps, cpus=cpus,
                     model='spde')
  nrepars <- length(obj$env$parList()$x)
  stats <- rbind(stats, cbind(model='spde', ndata=ndata,
                              nrepars=nrepars, get_stats(fits)))
  saveRDS(stats, file='results/spde_stats.RDS')
}
