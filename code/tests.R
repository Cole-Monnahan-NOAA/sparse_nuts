


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
