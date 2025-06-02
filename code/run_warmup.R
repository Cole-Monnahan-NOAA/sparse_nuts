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
warmups.irt_2pl <- fit_warmups(obj.irt_2pl, model='irt_2pl', init='last.par.best')
warmups.irt_2pl_nc <- fit_warmups(obj.irt_2pl_nc, model='irt_2pl_nc', init='last.par.best')
