source(here::here('code/startup.R'))
source('code/load_TMB_objects.R')

## Warmup tests. Run 10 independent chains with warmup=1000 for
## 'auto' metric and see how long it takes to stabilize. I do 10
## chains b/c inits is not unique when using chains>1 due to
## limitation of stan_sample.
reps <- 1:10
cpus <- 10
warmups.sdmTMB <- fit_warmups(obj.sdmTMB, cpus=cpus, reps=reps)
warmups.pollock <- fit_warmups(obj.pollock, cpus=cpus, reps=reps)
warmups.wildf <- fit_warmups(obj.wildf, cpus=cpus, reps=reps)
warmups.swallows <- fit_warmups(obj.swallows, cpus=cpus, reps=reps)
warmups.dlm <- fit_warmups(obj.dlm, model='dlm', cpus=cpus, reps=reps)
warmups.wham <- fit_warmups(obj.wham, cpus=cpus, reps=reps)
warmups.sam <- fit_warmups(obj.sam, cpus=cpus, reps=reps)
warmups.petrel <- fit_warmups(obj.petrel, cpus=cpus, reps=reps)
warmups.gp_pois_regr <- fit_warmups(obj.gp_pois_regr, cpus=cpus, reps=reps)
warmups.salamanders <- fit_warmups(obj.salamanders, cpus=cpus, reps=reps)

source('code/load_RTMB_objects.R')
## Warmup tests. Run 10 independent chains with warmup=1000 for
## 'auto' metric and see how long it takes to stabilize. I do 10
## chains b/c inits is not unique when using chains>1 due to
## limitation of stan_sample.
reps <- 1:10
cpus <- 10
warmups.schools <- fit_warmups(obj.schools, model='schools', cpus=cpus, reps=reps)
warmups.diamonds <- fit_warmups(obj.diamonds, model='diamonds', cpus=cpus, reps=reps)
warmups.radon <- fit_warmups(obj.radon, model='radon', cpus=cpus, reps=reps)
warmups.kilpisjarvi <- fit_warmups(obj.kilpisjarvi, model='kilpisjarvi', cpus=cpus, reps=reps)
warmups.causal <- fit_warmups(obj.causal, model='causal', cpus=cpus, reps=reps)
warmups.irt_2pl <- fit_warmups(obj.irt_2pl, model='irt_2pl', cpus=cpus, reps=reps)
warmups.irt_2pl_nc <- fit_warmups(obj.irt_2pl_nc, model='irt_2pl_nc', cpus=cpus, reps=reps)
warmups.lynx <- fit_warmups(obj.radon, model='lynx', cpus=cpus, reps=reps)
