source('code/load_TMB_objects.R')
reps <- 1:3 # vector of replicate analyses
cpus <- 1 # parallel sessions for parallel chains!


# fits.simple <- fit_models(obj.simple, iter=2000,
#                           cpus=cpus, replicates=reps)
fits.sdmTMB <- fit_models(obj.sdmTMB, iter=2000,
                          control=list(adapt_delta=.9),
                          cpus=cpus, replicates=reps)
fits.pollock <- fit_models(obj.pollock, iter=2000,
                           init='last.par.best',
                           cpus=cpus, replicates=reps)
fits.wildf <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         init='last.par.best',
                         replicates=reps, control=list(adapt_delta=.95))
# turn on metric adaptation
fits.wildf_adapted <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         replicates=reps, control=list(adapt_delta=.95),
                         init='last.par.best',
                         adapt_metric = TRUE, warmup=1000, model='wildf_adapted')
# this one has some inits fail so fall back to mode
fits.swallows <- fit_models(obj.swallows, iter=2000, cpus=cpus,
                            replicates=reps, init='last.par.best',
                            control=list(adapt_delta=.99))
fits.dlm <- fit_models(obj.dlm, iter=2000, cpus=cpus,
                       replicates=reps, model='dlm')
fits.salamanders <- fit_models(obj.salamanders, iter=2000, cpus=cpus,
                       replicates=reps, model='salamanders')
fits.gp_pois_regr <- fit_models(obj.gp_pois_regr, iter=2000, cpus=cpus,
                               replicates=reps, model='gp_pois_regr',
                               control=list(adapt_delta=.99))
fits.petrel <- fit_models(obj.petrel, iter=2000, cpus=cpus,
                          init='last.par.best',
                          replicates=reps, model='petrel')
fits.wham <- fit_models(obj.wham, iter=2000, cpus=cpus,
                        init='last.par.best',
                        control=list(adapt_delta=.9),
                        replicates=reps, model='wham')
fits.sam <- fit_models(obj.sam, iter=2000, cpus=cpus,
                        init='last.par.best',
                        control=list(adapt_delta=.9),
                        replicates=reps, model='sam')

# Rerun with laplace turned on, called "embedded laplace
# approximation" (ELA) by Margossian et al
fits.sdmTMB_ELA <- fit_models(obj.sdmTMB, iter=2000,
                          control=list(adapt_delta=.9),
                          cpus=cpus, replicates=reps,
                          laplace=TRUE,
                          model='sdmTMB_ELA')
# fits.pollock_ELA <- fit_models(obj.pollock, iter=2000,
#                            cpus=cpus, replicates=reps,
#                            model='pollock_ELA', laplace=TRUE)
fits.wildf_ELA <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         init='last.par.best',
                         replicates=reps, control=list(adapt_delta=.95),
                         model='wildf_ELA', laplace=TRUE)
# this one has some inits fail so fall back to mode
fits.swallows_ELA <- fit_models(obj.swallows, iter=2000, cpus=cpus,
                            replicates=reps, init='last.par.best',
                            control=list(adapt_delta=.9),
                            laplace=TRUE, model='swallows_ELA')
fits.dlm_ELA <- fit_models(obj.dlm, iter=2000, cpus=cpus,
                       laplace=TRUE,
                       replicates=reps, model='dlm_ELA')
fits.salamanders_ELA <- fit_models(obj.salamanders, iter=2000, cpus=cpus,
                               replicates=reps, model='salamanders_ELA',
                               laplace=TRUE)
fits.gp_pois_regr_ELA <- fit_models(obj.gp_pois_regr, iter=2000, cpus=cpus,
                                warmup=1000,
                                replicates=reps, model='gp_pois_regr_ELA',
                                control=list(adapt_delta=.99),
                                laplace=TRUE)
# why is this sooo slow??
# 1000 transitions using 10 leapfrog steps per transition would take 10853.6 seconds.
# fits.petrel_ELA <- fit_models(obj.petrel, iter=2000, cpus=cpus,
#                                laplace=TRUE, init='last.par.best',
#                            replicates=reps, model='petrel_ELA')
## crazy slow again
# fits.wham_ELA <- fit_models(obj.wham, iter=2000, cpus=cpus,
#                         init='last.par.best', laplace=TRUE,
#                         control=list(adapt_delta=.9),
#                         replicates=reps, model='wham_ELA')
fits.sam_ELA <- fit_models(obj.sam, iter=2000, cpus=cpus,
                       init='last.par.best', laplace=TRUE,
                       control=list(adapt_delta=.9),
                       replicates=reps, model='sam_ELA')





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
fits.irt_2pl <- fit_models(obj.irt_2pl, iter=2000, cpus=cpus, replicates=reps,
                          model='irt_2pl', init='last.par.best',
                          globals=list(irt_2pl_dat=irt_2pl_dat),
                          control=list(adapt_delta=.9))
fits.irt_2pl_nc <- fit_models(obj.irt_2pl_nc, iter=2000, cpus=cpus, replicates=reps,
                           model='irt_2pl_nc', init='last.par.best',
                           globals=list(irt_2pl_dat=irt_2pl_dat),
                           control=list(adapt_delta=.9))
fits.lynx <- fit_models(obj.lynx, replicates=reps, cpus=cpus,
                         iter=2000, globals=globals.lynx,
                         control=list(adapt_delta=.9), init='last.par.best',
                         model='lynx')

# Rerun with laplace turned on, called "embedded laplace
# approximation" (ELA) by Margossian et al
fits.schools_ELA <- fit_models(obj.schools, iter=2000,
                           cpus=cpus, replicates=reps,
                           laplace=TRUE,
                           globals=list(schools_dat=schools_dat),
                           model='schools_ELA')
fits.diamonds_ELA <- fit_models(obj.diamonds, replicates=reps,
                            cpus=cpus, iter=2000,
                            laplace=TRUE, init='last.par.best',
                            globals=list(diamonds_dat=diamonds_dat),
                            model='diamonds_ELA')
fits.radon_ELA <- fit_models(obj.radon, replicates=reps, cpus=cpus,
                         iter=2000, globals=list(radon_dat=radon_dat),
                         control=list(adapt_delta=.9),
                         laplace=TRUE,
                         model='radon_ELA')
fits.causal_ELA <- fit_models(obj.causal, iter=2000, cpus=cpus, replicates=reps,
                          model='causal_ELA', init='last.par.best',
                          laplace=TRUE,
                          control=list(adapt_delta=.9))
fits.irt_2pl_ELA <- fit_models(obj.irt_2pl, iter=2000, cpus=cpus, replicates=reps,
                           model='irt_2pl_ELA', init='random',
                           laplace=TRUE,
                           globals=list(irt_2pl_dat=irt_2pl_dat),
                           control=list(adapt_delta=.9))
fits.irt_2pl_nc_ELA <- fit_models(obj.irt_2pl_nc, iter=2000, cpus=cpus, replicates=reps,
                               model='irt_2pl_nc_ELA', init='random',
                               laplace=TRUE,
                               globals=list(irt_2pl_dat=irt_2pl_dat),
                               control=list(adapt_delta=.9))
fits.lynx_ELA <- fit_models(obj.lynx, replicates=reps, cpus=cpus,
                            laplace=TRUE,
                            iter=2000, globals=globals.lynx,
                            control=list(adapt_delta=.9), init='last.par.best',
                            model='lynx_ELA')
