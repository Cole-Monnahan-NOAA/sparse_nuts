library(RTMB)

stats.sd <- list()
for(sd in 10^(0:8)){
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
  fits <- fit_models(obj, iter=2000, warmup=1000, cores=4,
                    chains=4, Q=Q,Qinv=M,
                    metric=c('unit', 'diag', 'dense'),
                    model='ratios', cpus=1,
                    globals=list(data=data), replicates=reps,
                    control=list(max_treedepth=15))
  stats.sd <- rbind(stats.sd, cbind(sd=sd, get_stats(fits)))
}
saveRDS(stats.sd, "results/ratio_stats.RDS")
saveRDS(fits, 'results/ratio_fits.RDS')



stats.cor <- list()
## want them piled up near 1
cors <- c(.01,1/(1+exp(-seq(-1,10, len=6))))
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
  fits <- fit_models(obj, iter=2000, warmup=1000, cores=4,
                    chains=4, Q=Q,Qinv=M, skip_optimization=TRUE,
                    metric=c('unit', 'diag', 'dense'),
                    model='cors', plot = FALSE,
                    globals=list(data=data), replicates=reps, cpus=1,
                    control=list(max_treedepth=15))
  stats.cor <- rbind(stats.cor, cbind(cor=cor, get_stats(fits)))
}
saveRDS(stats.cor, "results/cor_stats.RDS")
saveRDS(fits, 'results/cor_fits.RDS')
