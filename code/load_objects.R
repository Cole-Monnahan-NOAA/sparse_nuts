
schools_dat <- list(J = 8,
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
pars <- list(mu=0, logtau=0, eta=rep(1,8))
f <- function(pars){
  RTMB::getAll(schools_dat, pars)
  theta <- mu + exp(logtau) * eta;
  lp <- sum(dnorm(eta, 0,1, log=TRUE))+ # prior
         sum(dnorm(y,theta,sigma,log=TRUE))+ #likelihood
         logtau                          # jacobian
  return(-lp)
}
obj.schools <- RTMB::MakeADFun(func=f, parameters=pars,
                    random="eta", silent=TRUE)


TMB::runExample('simple')
obj.simple <- obj



library(dsem)
## library(MARSS)
## rm(list=ls())
## knitr::knit('C:/Users/cole.monnahan/dsem/vignettes/dynamic_factor_analysis.Rmd')
## obj <- mydsem_full$obj
## obj <- mydfa$obj
## saveRDS(obj, file='models/dsem_obj.RDS')

obj.dsem <- readRDS('models/dsem_obj.RDS')
obj.dsem$retape()

pk <- readRDS('models/pollock/fit.RDS')
compile('models/pollock/goa_pk_tmb.cpp')
dyn.load('models/pollock/goa_pk_tmb.dll')
obj.pollock <- MakeADFun(data=pk$input$dat, parameters=pk$input$pars,
                         random=pk$input$random, DLL='goa_pk_tmb',
                         map=pk$input$map, silent=TRUE)


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


dat <- readRDS(file='models/diamonds/dat.RDS')
Kc <- dat$K - 1;
Xc <- matrix(NA, nrow=dat$N, ncol=Kc) # centered version of X without an intercept
means_X <- rep(NA, Kc) # column means of X before centering
for(i in 2 : dat$K) {
  means_X[i - 1] = mean(dat$X[,i]);
  Xc[,i - 1] = dat$X[,i] - means_X[i - 1];
}
dat$means_X <- means_X
dat$Xc <- Xc
pars <- list(b=rep(0, Kc),## population-level effects
     Intercept=1,   ## temporary intercept for centered predictors
     logsigma=0) ## residual SD
f <- function(pars){
  getAll(pars,dat)
  sigma <- exp(logsigma)
  lp <- sum(dnorm(b, 0, 1, TRUE))
  lp <- lp+sum(dt((Intercept-8)/10, 3,log=TRUE) )
  ## target += student_t_lpdf(Intercept | 3, 8, 10);
  ##target += student_t_lpdf(sigma | 3, 0, 10) - 1 * student_t_lccdf(0 | 3, 0, 10);
  lp <- lp+sum(dt(sigma/10,df=3, log=TRUE)) #-    sum(pt(sigma/10,df=3, log=TRUE))
  ## likelihood including all constants
  if (!prior_only) {
    lp <- lp+sum(dnorm(Y, mean=Intercept + b*Xc, sd=sigma, log=TRUE))
    ##target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  }
  lp <- lp+logsigma ## Jacobian
  b_Intercept <- Intercept - sum(means_X*b)
  REPORT(b_Intercept)
  return(-lp)
}
obj.diamonds <- RTMB::MakeADFun(f, pars, random='b', silent=TRUE)
