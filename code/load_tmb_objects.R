#detach(package:RTMB)
library(TMB)


TMB::runExample('simple')
obj.simple <- obj



library(dsem)
# library(MARSS)
# rm(list=ls())
# knitr::knit('C:/Users/cole.monnahan/dsem/vignettes/dynamic_factor_analysis.Rmd')
# obj <- mydsem_full$obj
# obj <- mydfa$obj
# saveRDS(obj, file='models/dsem/dsem_obj.RDS')
obj.dsem <- readRDS('models/dsem/dsem_obj.RDS')
obj.dsem$retape()

pk <- readRDS('models/pollock/fit.RDS')
compile('models/pollock/pollock.cpp')
dyn.load('models/pollock/pollock.dll')
#pk$input$map$sigmaR <- factor(1)
obj.pollock <- MakeADFun(data=pk$input$dat, parameters=pk$input$pars,
                         random=c('dev_log_recruit', pk$input$random)[2],
                         DLL='pollock',
                         map=pk$input$map, silent=TRUE)
obj.pollock$par <- pk$opt$par
# pk$obj$retape()
# obj.pollock <- pk$obj
# opt <- with(obj.pollock, nlminb(pk$opt$par, fn,gr))
# opt <- with(obj.pollock, nlminb(opt$par, fn,gr))
# xx <- GOApollock::fit_pk(pk$input)
# sdrep <- sdreport(obj.pollock, getJointPrecision = TRUE)
# fittmp <- sample_sparse_tmb(obj.pollock, iter=5000, warmup=100, metric='sparse',
#                              chains=4, cores=4, control=list(metric='unit_e'))
# pairs_admb(fittmp, pars=1:6, order='slow')

library(sdmTMB)
## Example model
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)
obj <- fit$tmb_obj
## refit with ln_kappa mapped off since inestimable as is
pars2 <- fit$parlist
map2 <- fit$tmb_map
map2$ln_kappa <- factor(c(NA,NA))
obj.sdmtmb <- TMB::MakeADFun(data=obj$env$data, parameters=pars2,
                       map=map2, random=fit$tmb_random, silent=TRUE, DLL=obj$env$DLL)

data <- readRDS('models/swallows/data.RDS')
inits <- function()
    list(sigmayearphi=runif(1, 0, 2),
         sigmaphi=runif(1,0,2),
         sigmap=runif(1,0,2),
         a=rnorm(data$K-1, 3, 1), a1=rnorm(1, .5, 1),
         b0=rnorm(4, 3, sd=1), b1=rnorm(4, 0, .15),
         fameffphi_raw=rnorm(data$nfam,0,1),
         fameffp_raw=rnorm(data$nfam,0,1),
         yeareffphi_raw=rnorm(4, 0,1))
set.seed(23254)
inits <- inits()
compile('models/swallows/swallows.cpp')
dyn.load('models/swallows/swallows')
obj.swallows <- MakeADFun(data=data, parameters=inits,
                 random=names(inits)[8:10],
                 DLL='swallows')

data <- readRDS('models/wildf/data.RDS')
inits <- function()
  list(yearInterceptSD = runif(1, .05, .15),
       plantInterceptSD = runif(1, .05, .15),
       plantSlopeSD = runif(1, .05, .15),
       intercept = rnorm(data$Nstage, 0, 1),
       slope = rnorm(1, 0, 1),
       yearInterceptEffect_raw= rnorm(data$Nyear, 0, 1),
       plantInterceptEffect_raw= rnorm(data$Nplant, 0, 1),
       plantSlopeEffect_raw= rnorm(data$Nplant, 0, 1))
set.seed(135213)
inits <- inits()
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj.wildf <- MakeADFun(data=data, parameters=inits,
                       random=names(inits)[6:8],
                       DLL='wildf')
