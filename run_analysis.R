## Main script to run the analysis.

## Started 6/2024 by Cole

## devtools::load_all('C:/Users/cole.monnahan/adnuts/')
## devtools::install('C:/Users/cole.monnahan/adnuts/')
## devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
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
cpus <- 3
reps <- 1:6
source("code/run_timings.R")

## increasing dimensionality of a SPDE object
cpus <- 3
reps <- 1:3
source('code/run_spde.R')

detach(package:TMB)
detach(package:RTMB)

source('code/load_TMB_objects.R')
cpus <- 1; reps <- 1:3
fits.simple <- fit_models(obj.simple, iter=2000, cpus=cpus, replicates=reps,
                          metrics=c('auto', 'unit'))
fits.sdmtmb <- fit_models(obj.sdmtmb, iter=2000, cpus=cpus, replicates=reps,
                          metrics=c('auto', 'unit'))
fits.pollock <- fit_models(obj.pollock, iter=2000, cpus=cpus, replicates=reps,
                           metrics=c('auto', 'unit'))
fits.wildf <- fit_models(obj.wildf, iter=2000, cpus=cpus,
                         replicates=reps, control=list(adapt_delta=.95),
                         metrics=c('auto', 'unit'))
fits.swallows <- fit_models(obj.swallows, iter=2000, cpus=cpus, replicates=reps,
                            metrics=c('auto', 'unit'),
                            control=list(adapt_delta=.99))
fits.dsem <- fit_models(obj.dsem, iter=2000, cpus=cpus, replicates=reps,
                        metrics=c('auto', 'unit'),
                        control=list(adapt_delta=.95))

detach(package:TMB)
source('code/load_RTMB_objects.R')
reps <- 1
cpus <- 1
fits.schools <- fit_models(obj.schools, iter=2000, cpus=cpus, replicates=reps,
                           globals=list(schools_dat=schools_dat),
                           model='schools')
fits.diamonds <- fit_models(obj.diamonds, cpus=cpus, iter=2000,
                            globals=list(diamonds_dat=diamonds_dat),
                            model='diamonds')
fits.gp_pois <- fit_models(obj.gp_pois, cpus=cpus, iter=2000,
                           globals=list(diamonds_dat=diamonds_dat),
                           control=list(adapt_delta=.95),
                           model='gp_pois')
fits.kilpisjarvi <- fit_models(obj.kilpisjarvi, cpus=cpus, iter=2000,
                           globals=list(kilpisjarvi_dat=kilpisjarvi_dat),
                           control=list(adapt_delta=.9),
                           metrics = c('unit', 'diag', 'dense'),
                           model='kilpisjarvi')




ff <- list.files('results',  pattern='_fits.RDS', full.names = TRUE)
results <- ff |> sapply(readRDS)
stats <- lapply(results, get_stats) |> bind_rows() |>
  group_by(model) |>
  mutate(rel_eff=eff/mean(eff[metric=='unit'])) |> ungroup()
stats2 <- stats |> filter(metric %in% c('sparse', 'dense')) |>
  filter( !(metric=='sparse' & model=='diamonds'))
g <- ggplot(stats2, aes(x=model, y=rel_eff)) + geom_jitter(width=.1, height=0) +
  scale_y_log10() + geom_hline(yintercept=1) +
  labs(x=NULL, y='Efficiency (minESS/time) relative to unit metric')
g
ggsave('plots/relative_perf.png', g, width=6, height=4)

fits <- readRDS('results/swallows_fits.RDS')
mean(fits[[1]]$mle$Q==0)
tmp <- fits[[1]]$mle$Q
tmp[tmp!=0] <- 1
Matrix::image(tmp)
tmp <- as.matrix(tmp)

mean(tmp==0)
test <- sample_sparse_tmb(obj.swallows, iter=5000, warmup=1000,
                          metric='dense', chains=6, cores=6,
                          control=list(adapt_delta=.99))
pars <- c('sigmayearphi', 'yeareffphi_raw[1]','yeareffphi_raw[2]','yeareffphi_raw[3]','yeareffphi_raw[4]')[1:3]
pairs_admb(fits[[3]],  pars=pars)

fits <- readRDS('results/wildf_fits.RDS')
pars <- c('slope', 'plantSlopeSD', 'plantSlopeEffect_raw[457]')
pairs_admb(fits[[1]],  pars=pars)
pairs_admb(fits[[2]],  pars=pars)
pairs_admb(fits[[1]], pars=1:8, order='mismatch')




## validate the model by running it on known distributions and
## checking for the correct properties
library(RTMB)
set.seed(23)
cov <- rWishart(n=1, df=2, Sigma=diag(2))[,,1]
cov2cor(cov)
pars <- list(x=c(0,0))
dat <- list(Sigma=cov)
f <- function(pars){
  getAll(pars,dat)
  -dmvnorm(x=x, mu=c(1,1), Sigma = Sigma, log=TRUE)
}
f <- function(pars){
  getAll(pars)
  -sum(dnorm(x=x, mean=c(1,1), sd = c(1,1), log=TRUE))
}
obj <- MakeADFun(f, pars, silent=TRUE)
obj$fn()
obj$gr()
opt <- with(obj, nlminb(par,gr))
sdrep <- sdreport(obj)
sdrep$cov.fixed

fits.validate <- fit_models(obj, iter=2000000, replicates=1, thin=10,
                            warmup=1000, cores=3, chains=3, globals=list(dat=dat),
                            metrics=c('unit', 'diag', 'dense')[1])
posts <- lapply(fits.validate, function(x) cbind(metric=x$metric, as.data.frame(x))) |>
  bind_rows() %>% setNames(c('metric', 'x1', 'x2')) %>%
  group_by(metric) %>%
  mutate(iter=1:n(), x1mean=cummean(x1)-1, x2mean=cummean(x2)-1) |> ungroup()

ggplot(filter(posts, iter>1000), aes(iter, abs(x1mean), color=metric)) + geom_line() + scale_y_log10() +
  scale_x_log10()

ggplot(filter(posts, iter>1000), aes(iter, x1mean, color=metric)) + geom_line() +
  scale_x_log10()


group_by(posts, metric) %>% filter(iter==max(iter))




library(posteriordb)
library(tidyr); library(dplyr); library(ggplot2)
#my_pdb <- pdb_local(path='C:/Users/cole.monnahan/posteriordb-r/')
my_pdb <- pdb_github()
pnms <- posterior_names(my_pdb)
po <- posterior(pnms[i], my_pdb)

stan_code('mcycle_gp-accel_gp')
stan_data('hmm_example')

out <- lapply(pnms, function(nm){
  ref <- tryCatch(reference_posterior_draws(nm), error=function(e) 'none')
  if(is.character(ref)) return(NULL)
  df <- lapply(ref[1:length(ref)], \(x) as.data.frame(x)) |> bind_rows()
  cor <- cor(df)
  data.frame(name=nm, npars=nrow(cor), cor=cor[upper.tri(cor, diag=FALSE)])
})
out <- do.call(rbind, out)
ggplot(out, aes(name, y=abs(cor))) + geom_violin() + coord_flip()

outsum <- out |> group_by(name) |> summarize(npars=npars[1], maxcor=max(abs(cor)))
outsum |> arrange(desc(maxcor))
outsum |> arrange(desc(npars))
ggplot(outsum, aes(npars, maxcor)) + geom_point()
