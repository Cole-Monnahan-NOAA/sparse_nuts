
library(glmmTMB)
data(Salamanders)
dat <- Salamanders

metrics <- c('unit', 'diag', 'dense', 'sparse')
stats <- list()
for(nreps in c(1,4,16,64,256, 512)){
  set.seed(nreps)
  message("Starting analysis for nreps=",nreps)
  repid <- rep(1:nreps, times=ndat)
  siteid <- rep(1:ndat, each=nreps)
  dat2 <- dat[siteid,]
  dat2$site2 <- as.factor(paste0(dat2$site, '-',repid))
  # checks that it worked
  stopifnot(nrow(dat2)/nrow(dat)==nreps)
  stopifnot(length(levels(dat2$site2))/length(levels(dat$site))==nreps)
  obj <- glmmTMB(count ~ spp * mined + (1 | site2), data = dat2, family="nbinom2")$obj
  fits <- fit_models(obj, chains=4, cores=4,
                     metrics=metrics,
                     iter=2000,
                     init='random', replicates=reps, cpus=cpus,
                     model='glmmTMB', plot=FALSE)
  nrepars <- length(obj$env$random)
  stats <- rbind(stats, cbind(nreps=nreps, nrepars=nrepars, get_stats(fits)))
  saveRDS(stats, file='results/glmmTMB_stats.RDS')
  g <- ggplot(stats, aes(nrepars, eff, group=interaction(metric,replicate), color=metric)) +
    geom_point(alpha=.5) +
    geom_line() +
    scale_y_log10() + labs(x='# of random effects', y='Efficiency (ESS/t)')
  ggsave('plots/glmmTMB_stats.png', g, width=7, height=5)
  print(g)
}

