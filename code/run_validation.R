

# Check that all metrics give the same posteriors
library(RTMB)
library(adnuts)
setwd(here::here())
source(here::here('code/startup.R'))

metrics <- c('unit', 'diag', 'dense', 'sparse', 'RTMBtape',  'sparse-J')

source('code/load_rtmb_objects.R')
# First on a known simple example
pars <- list(xnorm=1,logxgamma=1, xmvn=c(1,1))
func <- function(pars){
  getAll(pars)
  lp <-  logxgamma +
    dnorm(x=xnorm, mean=1.11, sd=.5, log=TRUE) +
    dgamma(x=exp(logxgamma), shape=.98, scale=.1, log=TRUE) +
    RTMB::dmvnorm(x=xmvn, mu=c(0,1), Sigma = matrix(c(1,.9,.9,1), nrow=2), log=TRUE)
  return(-lp)
}
func(pars)
obj <- MakeADFun(func, parameters = pars, random='xmvn')
obj$fn()
opt <- with(obj, nlminb(par, fn, gr))
thin <- 10
iter <- 4100*thin
wm <- thin*100
fits <- fit_models(obj, warmup=wm, iter=iter, replicates=1,
                   chains=1, cores=1, init='random',
                   cpus=1, thin=thin, metrics = metrics, plot = FALSE)
lapply(fits, \(x) x)
posts <- lapply(fits, \(x)
                cbind(metric=x$metric,
                      extract_samples(x, inc_lp=TRUE, inc_warmup = FALSE))) |>
  bind_rows()
group_by(posts, metric) |> slice_head(n=1)
n <- fits[[1]]$iter * dim(fits[[1]]$samples)[2]
set.seed(12314)
tmp <- mvtnorm::rmvnorm(n, mean=c(0,1), sigma = matrix(c(1,.9,.9,1), nrow=2))
postiid <- data.frame(metric='iid', xnorm=rnorm(n, 1.11, .5),
                      logxgamma=log(rgamma(n, shape=.98, scale=.1)),
                      'xmvn[1]'=tmp[,1], 'xmvn[2]'=tmp[,2], lp__=1) |>
  setNames(names(posts))
postiid$lp__ <-
  RTMB::dmvnorm(tmp, mu=c(0,1), Sigma = matrix(c(1,.9,.9,1), nrow=2), log=TRUE)+
  dnorm(postiid$xnorm, 1.11,.5, log=TRUE)+
  dgamma(exp(postiid$logxgamma), shape=.98, scale=.1, log=TRUE) +
  postiid$logxgamma
#postiid$lp__ <- NA
posts <- bind_rows(posts, postiid)
pl <- pivot_longer(posts, -c(metric)) |>
  group_by(metric, name) |>
  mutate(iter=1:n(), cummean=cummean(value)) |> ungroup()

g <- ggplot(filter(pl, iter>50), aes(iter, cummean, color=metric)) +
  facet_wrap('name', scales='free', nrow=3) + geom_line() +
  labs(y='Cumulative mean', x='Iteration')
ggsave('plots/validation_known.png', g, width=7, height=5)


if('RTMB' %in% .packages()) detach(package:RTMB)
library(TMB)
library(adnuts)
TMB::runExample('simple')
thin <- 10
iter <- 2100*thin
wm <- thin*100
fits <- fit_models(obj, warmup=wm, iter=iter, replicates=1,
                   cpus=1, chains=4, cores=1, init='last.par.best', plot=FALSE,
                   thin=thin, metrics = metrics, control=list(max_treedepth=6))
lapply(fits, \(x) x)
pars <- c('u[1]', 'u[114]', 'beta[1]', 'beta[2]', 'logsdu', 'logsd0', 'lp__')
posts <- lapply(fits, \(x)
                cbind(metric=x$metric,
                      extract_samples(x, inc_lp=TRUE, inc_warmup = FALSE))) |>
  bind_rows() |> select(all_of(c('metric',pars)))

pl <- pivot_longer(posts, -c(metric)) |>
  group_by(metric, name) |>
  mutate(iter=1:n(), cummean=cummean(value)) |> ungroup()

g <- ggplot(filter(pl, iter>1000), aes(iter, cummean, color=metric)) +
  facet_wrap('name', scales='free', nrow=3) + geom_line() +
  labs(y='Cumulative mean', x='Iteration')
ggsave('plots/validation_simple.png', g, width=7, height=5)
