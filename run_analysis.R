## Main script to run the analysis.

## Started 6/2024 by Cole

## devtools::load_all('C:/Users/cole.monnahan/adnuts/')
## devtools::install('C:/Users/cole.monnahan/adnuts/')
## devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
setwd(here::here())
source("code/startup.R")

# ## Demo with simple
# TMB::runExample('simple')
# ## TMB::runExample('ar1xar1')
# fit.demo <- sample_sparse_tmb(obj, iter=500, warmup=100, init='random',
#                               chains=2, cores=1, metric='sparse', seed=124)
# launch_shinyadmb(fit.demo)



## Run benchmark and timings on spatial Poisson model.
## source("code/run_spatial_benchmark.R")

## Efficiency across a bivariate normal with varying marginal SDs and correlations
cpus <- 3
reps <- 1:6
source("code/run_timings.R")

## increasing dimensionality of a SPDE object
cpus <- 3
reps <- 1:3
source('code/run_spde.R')

## increasing dimensionality of negative binomial regression
cpus <- 3
reps <- 1:3
source('code/run_glmmTMB.R')

# increasing dimensionality of vector autoregressive model in
# 'dsem', fitted to original wolf-moose model from
# https://james-thorson-noaa.github.io/dsem/articles/vignette.html#comparison-with-vector-autoregressive-models
# but modified to have a normal family w/ fixed observation
# SD=0.1; then simulated data by resimulating the GMRF and new
# data.
cpus <- 1
reps <- 1:3
source('code/run_VAR.R')

source('code/load_TMB_objects.R')

## Warmup tests. Run 10 independent chains with warmup=1000 for
## 'auto' metric and see how long it takes to stabilize. I do 10
## chains b/c inits is not unique when using chains>1 due to
## limitation of stan_sample.
reps <- 1:10
cpus <- 5
warmups.simple <- fit_warmups(obj.simple)
warmups.sdmTMB <- fit_warmups(obj.sdmTMB)
warmups.pollock <- fit_warmups(obj.pollock, init='last.par.best')
warmups.wildf <- fit_warmups(obj.wildf, init='last.par.best')
warmups.swallows <- fit_warmups(obj.swallows, init='last.par.best')
warmups.dlm <- fit_warmups(obj.dlm, model='dlm')


reps <- 1:3 # vector of replicate analyses
cpus <- 1 # parallel sessions for parallel chains!
fits.simple <- fit_models(obj.simple, iter=2000,
                          cpus=cpus, replicates=reps)
fits.sdmTMB <- fit_models(obj.sdmTMB, iter=2000,
                          control=list(adapt_delta=.9),
                          cpus=cpus, replicates=reps)
fits.pollock <- fit_models(obj.pollock, iter=2000,
                           cpus=cpus, replicates=reps)
fits.wildf <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         replicates=reps, control=list(adapt_delta=.95))
# this one has some inits fail so fall back to mode
fits.swallows <- fit_models(obj.swallows, iter=2000, cpus=cpus,
                            replicates=reps, init='last.par.best',
                            control=list(adapt_delta=.99))
fits.dlm <- fit_models(obj.dlm, iter=2000, cpus=cpus,
                       replicates=reps, model='dlm')


source('code/load_RTMB_objects.R')
## Warmup tests. Run 10 independent chains with warmup=1000 for
## 'auto' metric and see how long it takes to stabilize. I do 10
## chains b/c inits is not unique when using chains>1 due to
## limitation of stan_sample.
reps <- 1:10
cpus <- 5
warmups.schools <- fit_warmups(obj.schools, model='schools')
warmups.diamonds <- fit_warmups(obj.diamonds, model='diamonds')
warmups.radon <- fit_warmups(obj.radon, model='radon')
warmups.kilpisjarvi <- fit_warmups(obj.kilpisjarvi, model='kilpisjarvi')
warmups.causal <- fit_warmups(obj.causal, model='causal', init='last.par.best')


reps <- 1:3 # vector of replicate analyses
cpus <- 1 # parallel sessions for parallel chains!
fits.schools <- fit_models(obj.schools, iter=2000,
                           cpus=cpus, replicates=reps,
                           globals=list(schools_dat=schools_dat),
                           model='schools')
fits.diamonds <- fit_models(obj.diamonds, replicates=reps,
                            cpus=cpus, iter=2000,
                            globals=list(diamonds_dat=diamonds_dat),
                            model='diamonds')
# fits.gp_pois <- fit_models(obj.gp_pois, cpus=cpus,
#                            replicates=reps, iter=2000,
#                            globals=list(gp_pois_dat=gp_pois_dat),
#                            control=list(adapt_delta=.95),
#                            model='gp_pois')
fits.radon <- fit_models(obj.radon, replicates=reps, cpus=cpus,
                         iter=2000, globals=list(radon_dat=radon_dat),
                         control=list(adapt_delta=.9),
                         model='radon')
fits.kilpisjarvi <- fit_models(obj.kilpisjarvi, replicates=reps,
                               cpus=cpus, iter=2000,
                               globals=list(kilpisjarvi_dat=kilpisjarvi_dat),
                               control=list(adapt_delta=.9),
                               model='kilpisjarvi')
fits.causal <- fit_models(obj.causal, iter=2000, cpus=cpus, replicates=reps,
                          model='causal', init='last.par.best',
                          control=list(adapt_delta=.9))




results <- list.files('results',  pattern='_fits.RDS', full.names = TRUE) |> lapply(readRDS)
stats <- lapply(results, get_stats) |> bind_rows() |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit']),
         mean_rel_eff=mean(rel_eff)) |> ungroup() |>
  arrange(desc(mean_rel_eff)) |>
  mutate(model=factor(model, levels=unique(model)),
         metric2=metricf2(metric))
g <- ggplot(stats, aes(x=model, y=rel_eff, color=metric2)) +
  geom_jitter(width=.1, height=0) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_log10() + geom_hline(yintercept=1) +
    labs(x=NULL, y='Efficiency (minESS/time) relative to unit metric')
g
ggsave('plots/relative_perf.png', g, width=6, height=4)


warmups <- lapply(list.files('results/warmups/', full.names = TRUE), readRDS)
sp <- lapply(warmups, function(y){
  lapply(y, function(x){
  cat(x$model, x$metric, x$replicate, "\n")
  cbind(model=x$model, metric=x$metric, replicate=x$replicate, plot_sampler_params(x, plot=FALSE)$data )
  }
)
}
) |> bind_rows()
eps <- filter(sp, variable=='log_stepsize' & iteration < 1000)
g <- ggplot(eps, aes(iteration, y=value,
                     group=interaction(chain,replicate, metric))) +
  geom_line(alpha=.1) +
  facet_wrap('model', ncol=3,scale='free_y') + labs(y='log step-size')
ggsave('plots/warmup.png',g, width=8, height=8)




