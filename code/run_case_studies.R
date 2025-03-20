source('code/load_TMB_objects.R')
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
fits.salamanders <- fit_models(obj.salamanders, iter=2000, cpus=cpus,
                       replicates=reps, model='salamanders')



source('code/load_RTMB_objects.R')
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
