
dir.create('results/case_studies', showWarnings=FALSE)
dir.create('plots/case_studies', showWarnings=FALSE)
do.tmbstan <- TRUE
# save files to separate folder
out <- 'results/case_studies'
reps <- 1:3 # vector of replicate analyses
cpus <- 1 # parallel sessions for parallel chains!


source('code/load_TMB_objects.R')
fits.pollock <- fit_models(obj.pollock,
                           do.tmbstan=do.tmbstan,
                           cpus=cpus, replicates=reps,
                           outpath=out)
fits.wildf <- fit_models(obj.wildf,  cpus=cpus,
                         do.tmbstan=do.tmbstan,
                        # init='last.par.best',
                         replicates=reps, control=list(adapt_delta=.95),
                         outpath=out)
# turn on metric adaptation and longer warmup
fits.wildf_adapted <- fit_models(obj.wildf,  cpus=cpus,
                                 replicates=reps, control=list(adapt_delta=.95),
                                 adapt_metric=TRUE, num_warmup=1000,
                                 do.tmbstan=do.tmbstan,
                                 model='wildf_adapted',
                                 outpath=out)
# this one has some inits fail so fall back to mode
fits.swallows <- fit_models(obj.swallows,  cpus=cpus,
                            replicates=reps,
                            do.tmbstan=do.tmbstan,
                            #init='last.par.best',
                            control=list(adapt_delta=.99),
                            outpath=out)
fits.dlm <- fit_models(obj.dlm,  cpus=cpus,
                       replicates=reps, model='dlm',
                       do.tmbstan=do.tmbstan,
                       control=list(adapt_delta=.95),
                       outpath=out)
fits.salamanders <- fit_models(obj.salamanders,  cpus=cpus,
                               do.tmbstan=do.tmbstan,
                               control=list(adapt_delta=.9),
                               replicates=reps, model='salamanders',
                               outpath=out)
fits.gp_pois_regr <- fit_models(obj.gp_pois_regr,  cpus=cpus,
                               replicates=reps, model='gp_pois_regr',
                               do.tmbstan=do.tmbstan,
                               control=list(adapt_delta=.99),
                               outpath=out)
fits.petrel <- fit_models(obj.petrel,  cpus=cpus,
                          do.tmbstan=do.tmbstan,
                        #  init='last.par.best',
                          replicates=reps, model='petrel',
                          outpath=out)
fits.wham <- fit_models(obj.wham,  cpus=cpus,
                      #  init='last.par.best',
                      do.tmbstan=do.tmbstan,
                        control=list(adapt_delta=.9),
                        replicates=reps, model='wham',
                        outpath=out)
fits.sam <- fit_models(obj.sam,  cpus=cpus,
                     # throws error about tmb_forward https://github.com/kaskr/tmbstan/issues/11
                     do.tmbstan=FALSE,
                     control=list(adapt_delta=.9),
                     replicates=reps, model='sam',
                     outpath=out)
fits.sdmTMB <- fit_models(obj.sdmTMB,
                          control=list(adapt_delta=.9),
                          cpus=cpus, replicates=1,
                          do.tmbstan=FALSE,
                          outpath=out)
# Rerun with laplace turned on, called "embedded laplace
# approximation" (ELA) by Margossian et al

# fits.pollock_ELA <- fit_models(obj.pollock,
#                            cpus=cpus, replicates=reps,
#                            model='pollock_ELA', laplace=TRUE)
fits.wildf_ELA <- fit_models(obj.wildf,  cpus=cpus,
                         #init='last.par.best',
                         replicates=reps, control=list(adapt_delta=.95),
                         model='wildf_ELA', laplace=TRUE,
                         outpath=out)
# this one has some inits fail so fall back to mode
fits.swallows_ELA <- fit_models(obj.swallows,  cpus=cpus,
                            replicates=reps, init='last.par.best',
                            control=list(adapt_delta=.9),
                            laplace=TRUE, model='swallows_ELA',
                            outpath=out)
fits.dlm_ELA <- fit_models(obj.dlm,  cpus=cpus,
                       laplace=TRUE,
                       replicates=reps, model='dlm_ELA',
                       outpath=out)
fits.salamanders_ELA <- fit_models(obj.salamanders,  cpus=cpus,
                               replicates=reps, model='salamanders_ELA',
                               laplace=TRUE,
                               outpath=out)
fits.gp_pois_regr_ELA <- fit_models(obj.gp_pois_regr,  cpus=cpus,
                                warmup=1000,
                                replicates=reps, model='gp_pois_regr_ELA',
                                control=list(adapt_delta=.99),
                                laplace=TRUE, outpath=out)
# why is this sooo slow??
# 1000 transitions using 10 leapfrog steps per transition would take 10853.6 seconds.
# fits.petrel_ELA <- fit_models(obj.petrel,  cpus=cpus,
#                                laplace=TRUE, init='last.par.best',
#                            replicates=reps, model='petrel_ELA', outpath=out)
## crazy slow again
# fits.wham_ELA <- fit_models(obj.wham,  cpus=cpus,
#                         init='last.par.best', laplace=TRUE,
#                         control=list(adapt_delta=.9),
#                         replicates=reps, model='wham_ELA', outpath=out)
fits.sam_ELA <- fit_models(obj.sam,  cpus=cpus,
                       #init='last.par.best',
                       laplace=TRUE,
                       control=list(adapt_delta=.9),
                       replicates=reps, model='sam_ELA',
                       outpath=out)
fits.sdmTMB_ELA <- fit_models(obj.sdmTMB,
                              control=list(adapt_delta=.9),
                              cpus=cpus, replicates=reps,
                              laplace=TRUE,
                              model='sdmTMB_ELA',
                              outpath=out)


source('code/load_RTMB_objects.R')
reps <- 1:3 # vector of replicate analyses
cpus <- 1 # parallel sessions for parallel chains!
fits.schools <- fit_models(obj.schools,
                           cpus=cpus, replicates=reps,
                           globals=list(schools_dat=schools_dat),
                           control=list(adapt_delta=.95),
                           do.tmbstan=do.tmbstan,
                           model='schools',
                           outpath=out)
fits.diamonds <- fit_models(obj.diamonds, replicates=reps,
                            cpus=cpus,
                            do.tmbstan=do.tmbstan,
                            globals=list(diamonds_dat=diamonds_dat),
                            model='diamonds',
                            outpath=out)
fits.radon <- fit_models(obj.radon, replicates=reps, cpus=cpus,
                          globals=list(radon_dat=radon_dat),
                         control=list(adapt_delta=.9),
                         do.tmbstan=do.tmbstan,
                         model='radon',
                         outpath=out)
fits.kilpisjarvi <- fit_models(obj.kilpisjarvi, replicates=reps,
                               cpus=cpus,
                               do.tmbstan=do.tmbstan,
                               globals=list(kilpisjarvi_dat=kilpisjarvi_dat),
                               control=list(adapt_delta=.9),
                               model='kilpisjarvi',
                               outpath=out)
fits.causal <- fit_models(obj.causal,  cpus=cpus, replicates=reps,
                          model='causal',
                          # can't do tmbstan for this one since in package
                          do.tmbstan=FALSE,
                          #init='last.par.best',
                          control=list(adapt_delta=.9),
                          outpath=out)
fits.irt_2pl <- fit_models(obj.irt_2pl,  cpus=cpus, replicates=reps,
                          model='irt_2pl',
                          do.tmbstan=do.tmbstan,
                          #init='last.par.best',
                          globals=list(irt_2pl_dat=irt_2pl_dat),
                          control=list(adapt_delta=.9),
                          outpath=out)
fits.irt_2pl_nc <- fit_models(obj.irt_2pl_nc,  cpus=cpus, replicates=reps,
                           model='irt_2pl_nc',
                           do.tmbstan=do.tmbstan,
                           #init='last.par.best',
                           globals=list(irt_2pl_dat=irt_2pl_dat),
                           control=list(adapt_delta=.9),
                           outpath=out)
fits.lynx <- fit_models(obj.lynx, replicates=reps, cpus=cpus,
                        globals=globals.lynx,
                        do.tmbstan=do.tmbstan,
                        control=list(adapt_delta=.99),
                        #init='last.par.best',
                        model='lynx',
                        outpath=out)

# Rerun with laplace turned on, called "embedded laplace
# approximation" (ELA) by Margossian et al
fits.schools_ELA <- fit_models(obj.schools,
                           cpus=cpus, replicates=reps,
                           laplace=TRUE,
                           globals=list(schools_dat=schools_dat),
                           model='schools_ELA', outpath=out)
fits.diamonds_ELA <- fit_models(obj.diamonds, replicates=reps,
                            cpus=cpus,
                            laplace=TRUE, #init='last.par.best',
                            globals=list(diamonds_dat=diamonds_dat),
                            model='diamonds_ELA', outpath=out)
fits.radon_ELA <- fit_models(obj.radon, replicates=reps, cpus=cpus,
                          globals=list(radon_dat=radon_dat),
                         control=list(adapt_delta=.9),
                         laplace=TRUE,
                         model='radon_ELA',
                         outpath=out)
fits.causal_ELA <- fit_models(obj.causal,  cpus=cpus, replicates=reps,
                          model='causal_ELA', #init='last.par.best',
                          laplace=TRUE,
                          control=list(adapt_delta=.9),
                          outpath=out)
fits.irt_2pl_ELA <- fit_models(obj.irt_2pl,  cpus=cpus, replicates=reps,
                           model='irt_2pl_ELA',
                           laplace=TRUE,
                           globals=list(irt_2pl_dat=irt_2pl_dat),
                           control=list(adapt_delta=.9),
                           outpath=out)
fits.irt_2pl_nc_ELA <- fit_models(obj.irt_2pl_nc,  cpus=cpus, replicates=reps,
                               model='irt_2pl_nc_ELA',
                               laplace=TRUE,
                               globals=list(irt_2pl_dat=irt_2pl_dat),
                               control=list(adapt_delta=.9),
                               outpath=out)
fits.lynx_ELA <- fit_models(obj.lynx, replicates=reps, cpus=cpus,
                            laplace=TRUE,
                             globals=globals.lynx,
                            control=list(adapt_delta=.9), #init='last.par.best',
                            model='lynx_ELA', outpath=out)
