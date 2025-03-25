
library(glmmTMB)
data(Salamanders)
dat <- Salamanders
ndat <- nrow(dat)
nsites <- length(unique(dat$site))
ndat/nsites
mod = glmmTMB(count~spp * mined + (1|site), Salamanders, family="nbinom2")
########################################################
#Benchmark with data sets with more random effect levels


# code from Mollie Brooks from the glmmTMB paper
sims1 = lapply(simulate(mod, nsim=512, seed=111),
               function(count){
                 cbind(count, Salamanders[,c('site', 'mined', 'spp')])
               })

n = nrow(Salamanders)
nRE=length(unique(Salamanders$site))
bigdat0 = lapply(sitereps, function(x) do.call(rbind, sims1[1:x]))
bigdat =  lapply(1:length(sitereps), function(x)
  data.frame(bigdat0[[x]], "grp"=paste0(bigdat0[[x]]$site, rep(1:sitereps[x], each=n))))

options(future.globals.maxSize= 4000*1024^2) # 4000 MB limit
metrics <- c('unit', 'diag', 'dense', 'sparse')
stats <- list()
sitereps = c(2, 4,8, 16, 32,64, 128, 256)
for(ii in 5:length(sitereps)){
  nreps <- sitereps[ii]
  set.seed(nreps)
  message("Starting analysis for nreps=",nreps)
  dat2 <- bigdat[[ii]]
  obj <- glmmTMB(count ~ spp * mined + (1 | grp), data = dat2, family="nbinom2")$obj
  fits <- fit_models(obj, chains=4, cores=4,
                     metrics=metrics,
                     iter=2000,
                     init='random', replicates=reps,
                     cpus=ifelse(nreps<=128, cpus,1),
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
