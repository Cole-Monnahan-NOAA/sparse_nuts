## Main script to run the analysis.

## Started 6/2024 by Cole

## devtools::load_all('C:/Users/cole.monnahan/adnuts/')
## devtools::install('C:/Users/cole.monnahan/adnuts/')
library(adnuts)
source("code/startup.R")

TMB::runExample('simple')
## Demo with simple
## TMB::runExample('ar1xar1')
fit.demo <- sample_sparse_tmb(obj, iter=500, warmup=100, init='random',
                              chains=2, cores=1, metric='sparse', seed=124)
launch_shinyadmb(fit.demo)

## source("code/run_spatial_benchmark.R")

## Efficiency across a bivariate normal with varying marginal SDs and correlations
cpus <- 1
reps <- 1
source("code/run_timings.R")

## increasing dimensionality of a SPDE object
cpus <- 1
reps <- 1
source('code/run_spde.R')

source('code/load_objects.R')
cpus <- 1; reps <- 1:6
fits.schools <- fit_models(obj.schools, iter=2000, cpus=cpus, chains=1, replicates=reps)
fits.simple <- fit_models(obj.simple, iter=2000, cpus=cpus, replicates=reps)
fits.sdmtmb <- fit_models(obj.sdmtmb, iter=2000, cpus=cpus, replicates=reps)
fits.wildf <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         replicates=reps, control=list(adapt_delta=.9))
fits.swallows <- fit_models(obj.swallows, iter=2000, cpus=cpus, replicates=reps)
fits.pollock <- fit_models(obj.pollock, iter=2000, cpus=cpus, replicates=reps)
fits.diamonds <- fit_models(obj, cpus=1, iter=2000,
                            model='diamonds', globals=list(dat=dat))


