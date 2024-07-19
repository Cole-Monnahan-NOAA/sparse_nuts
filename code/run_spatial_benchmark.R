require(microbenchmark)
library(Matrix)
library(RTMB)
source("code/spatial_simulator.R")
## benchmark gradient calcs externally to show the impact of
## sparsity on the joint grad
bench.sparse <- bench.dense <- list()
ndata <- floor(sqrt(10^(1:3)))
ndata <- floor(sqrt(2^(4:11)))
for(n in ndata){
  print(n)
  obj <- sim_spde_dat(n, sparse=TRUE)
  nrepars <- length(obj$env$parList()$x)
  obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                          map=obj$env$map,
                          random=NULL, silent=TRUE,
                          DLL=obj$env$DLL)
  bench <- microbenchmark(obj2$gr(), times=50, unit='s')
  bench.sparse <- rbind(bench.sparse,
                        data.frame(model='sparse',
                                   ndata=n, nrepars=nrepars,
                                   time=summary(bench)$median[1]))
  obj <- sim_spde_dat(n, sparse=FALSE)
  obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                          map=obj$env$map,
                          random=NULL, silent=TRUE,
                          DLL=obj$env$DLL)
  bench <- microbenchmark(obj2$gr(), times=50, unit='s')
  bench.dense <- rbind(bench.dense,
                       data.frame(model='dense',
                                  ndata=n, nrepars=nrepars,
                                  time=summary(bench)$median[1]))
}
bench <- rbind(bench.sparse, bench.dense)
saveRDS(bench, 'results/bench_spde.RDS')

