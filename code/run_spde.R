options(future.globals.maxSize = 10000 * 1024^2) # 10 GB limit
stats <- list()
for(ndata in c(15,20,25,30,35,40,50,70,90,110)){
  set.seed(ndata)
  message("Starting analysis for ndata=",ndata)
  obj <- sim_spde_dat(ndata, TRUE, list(log_tau=factor(NA)))
  nrepars <- length(obj$env$last.par.best)
  fits <- fit_models(obj, iter=2000, chains=4, cores=4,
                     metrics = c('unit', 'auto'),
                     init='last.par.best',
                     replicates=reps,
                     cpus=ifelse(ndata>=40, 1, cpus),
                     model='spde', plot=FALSE)
  stats <- rbind(stats, cbind(ndata=ndata, nrepars=nrepars, get_stats(fits)))
  saveRDS(stats, file='results/spde_stats.RDS')
  g <- ggplot(stats, aes(nrepars, eff, group=interaction(metric,replicate), color=metric)) + geom_point(alpha=.5) +
    geom_line() +
    scale_y_log10() + labs(x='# of random effects', y='Efficiency (ESS/t)')
  ggsave('plots/spde_stats.png', g, width=7, height=5)
}

stats <- readRDS('results/spde_stats.RDS') %>% group_by(metric, nrepars) %>%
  mutate(mean_eff=mean(eff))
g <- ggplot(stats, aes(nrepars, mean_eff, color=metric, group=interaction(metric,replicate))) +
  geom_line() + scale_y_log10() + scale_x_log10()
g <- g + geom_jitter(mapping=aes(y=eff), width=.005)
ggsave('plots/spde_efficiency.png', g, width=7, height=5, units='in')
