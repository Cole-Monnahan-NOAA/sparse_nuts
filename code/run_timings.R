library(RTMB)

reps <- 1
cpus <- 1

stats.sd <- fits.list <- list()
for(sd in 10^(0:6)){
  Sigma <- matrix(c(1,0,0,sd), nrow=2)
  Mu <- rep(0,2)
  data <- list(Mu=Mu, Sigma=Sigma)
  params <- list(x=c(1,1))
  f <- function(params){
    ##  getAll(data,params)
    -sum(RTMB::dmvnorm(x=params$x, mu=data$Mu, Sigma=data$Sigma, log=TRUE))
  }
  obj <- RTMB::MakeADFun(func=f, parameters=params, silent=TRUE)
  opt <- with(obj, nlminb(par,fn,gr))
  sdr <- sdreport(obj)
  M <- sdr$cov.fixed
  Q <- solve(M)
  fits <- fit_models(obj, num_samples=1000, num_warmup=1000, cores=4,
                    chains=4, Q=Q,Qinv=M,
                    metric=c('stan', 'diag', 'dense'),
                    model='ratios', cpus=1, plot=FALSE,
                    init='unif',
                    globals=list(data=data), replicates=reps)
  fits.list <- c(fits.list, fits)
  stats.sd <- rbind(stats.sd, cbind(sd=sd, get_stats(fits)))
}
saveRDS(stats.sd, "results/ratio_stats.RDS")
saveRDS(fits.list, 'results/ratio_fits.RDS')



stats.cor <- fits.list <- list()
## want them piled up near 1
cors <- c(.01,1/(1+exp(-seq(-1,8, len=8))))
for(cor in cors){
  Sigma <- matrix(c(1,cor,cor,1), nrow=2)
  Mu <- rep(0,2)
  data <- list(Mu=Mu, Sigma=Sigma)
  params <- list(x=c(1,1))
  f <- function(params){
    ##  getAll(data,params)
    -sum(RTMB::dmvnorm(x=params$x, mu=data$Mu, Sigma=data$Sigma, log=TRUE))
  }
  obj <- RTMB::MakeADFun(func=f, parameters=params, silent=TRUE)
  opt <- with(obj, nlminb(par,fn,gr))
  sdr <- sdreport(obj)
  M <- sdr$cov.fixed
  Q <- solve(M)
  fits <- fit_models(obj, num_samples=1000, num_warmup=1000, cores=4,
                    chains=4, Q=Q,Qinv=M, skip_optimization=TRUE,
                    metric=c('stan', 'diag', 'dense'),
                    model='cors', plot = FALSE,
                    init='unif',
                    globals=list(data=data), replicates=reps, cpus=1)
  fits.list <- c(fits.list, fits)
  stats.cor <- rbind(stats.cor, cbind(cor=cor, get_stats(fits)))
}
saveRDS(stats.cor, "results/cor_stats.RDS")
saveRDS(fits.list, 'results/cor_fits.RDS')

lapply(fits.list, \(x) max(as.data.frame(extract_sampler_params(x))$treedepth)) |>
  as.numeric()
