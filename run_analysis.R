## Script to run the analysis for sparse NUTS manuscript (in prep).


## This will not work out of the box since many models need to be recompiled and
## linked locally.

## Started 6/2024 by Cole

## devtools::load_all('C:/Users/cole.monnahan/adnuts/')
## devtools::install('C:/Users/cole.monnahan/adnuts/')
## install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
## devtools::install_github('Cole-Monnahan-NOAA/SparseNUTS', ref='main')
setwd(here::here())
source("code/startup.R")

## Run benchmark and timings on spatial Poisson model.
## source("code/run_spatial_benchmark.R")

## Efficiency across a bivariate normal with varying marginal SDs and correlations
source("code/run_timings.R")

## SPDE benchmarking of gradient calls and memory usage
source('code/run_benchmarks.R')

## increasing dimensionality of a SPDE object
source('code/run_spde.R')

## increasing dimensionality of negative binomial regression
cpus <- 1
reps <- 1:3
source('code/run_glmmTMB.R')

# Run long warmups with adaptation of the Stan metric turned on to gauge if its
# needed
source("code/run_warmup.R")

# Run TMB and RTMB case studies
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



source("code/process_results.R")


# Need posteriors for this so do it last
source("code/run_pathfinder.R")



