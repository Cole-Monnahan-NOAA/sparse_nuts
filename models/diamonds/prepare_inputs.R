## remotes::install_github("stan-dev/posteriordb-r")
library(posteriordb)
my_pdb <- pdb_github()
po <- 'diamonds'
dat <- pdb_data(po)
code <- stan_code(po)
code
saveRDS(dat, file='models/diamonds/dat.RDS')


opt <- with(obj, nlminb(par,fn,gr))
library(adnuts)
library(StanEstimators)
fit <- sample_sparse_tmb(obj, iter=1000, warmup=200, chains=3,
                         cores=3, globals=list(dat=dat), metric='unit')

pairs_admb(fit, pars=1:10, order='slow')
plot_uncertainties(fit)
plot_sampler_params(fit)
