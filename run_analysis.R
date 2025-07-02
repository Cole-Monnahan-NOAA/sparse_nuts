## Main script to run the analysis.

## Started 6/2024 by Cole

## devtools::load_all('C:/Users/cole.monnahan/adnuts/')
## devtools::install('C:/Users/cole.monnahan/adnuts/')
## devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
setwd(here::here())
source("code/startup.R")

# ## Demo with simple
# library(TMB)
# library(adnuts)
# TMB::runExample('simple')
# fit <- sample_sparse_tmb(obj)
# launch_shinyadmb(fit)
# pairs_admb(fit, pars=1:5, order='slow')
# plot_marginals(fit, pars=1:5)
# plot_sampler_params(fit)
# plot_uncertainties(fit)


## Run benchmark and timings on spatial Poisson model.
## source("code/run_spatial_benchmark.R")

## Efficiency across a bivariate normal with varying marginal SDs and correlations
cpus <- 3
reps <- 1:3
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

source("code/run_warmup.R")
source('code/run_case_studies.R')


# Run some longer wildf ones to see if can identify why it's not
# working better
source("code/load_tmb_objects.R")
thin <- 3
warmup <- 2000
iter <- 3000*thin+warmup
fits.wildf2 <- fit_models(obj.wildf, warmup=warmup, replicates=1, cpus=1,
                          thin=thin, iter=iter, metrics=c('unit', 'auto'),
                          model = 'wildf2', init='last.par.best',
                          control=list(adapt_delta=.8))
fits.wildf3 <- fit_models(obj.wildf, warmup=warmup, replicates=1, cpus=1,
                          thin=thin, iter=iter, metrics=c('unit', 'auto'),
                          init = 'last.par.best',
                          control=list(adapt_delta=.999),
                          model = 'wildf3')




# Need posteriors for this so do it last
source("code/run_pathfinder.R")




