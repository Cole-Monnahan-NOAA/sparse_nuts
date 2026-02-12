setwd(here::here('models/wildf/'))

data <- readRDS('data.RDS')
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
compile('wildf.cpp')
dyn.load('wildf')
obj <- MakeADFun(data=data, parameters=inits,
                       random=names(inits)[6:8],
                       DLL='wildf')
saveRDS(obj, file='obj.wildf.RDS')

library(SparseNUTS)
fit <- sample_snuts(obj, iter=3000, warmup=500,
                         seed=1,
                         control=list(adapt_delta=.95))
pairs(fit, pars=1:5, order='slow')

obj$par |> length()
obj$env$par |> length() - obj$par |> length()
opt <- with(obj, nlminb(par,fn,gr))
Q <- sdreport(obj, getJointPrecision=TRUE)$jointPrecision
M <- solve(Q)
max(abs(cov2cor(M)[lower.tri(M)]))

post <- as.data.frame(fit)
cors <- cor(post)
max(abs(cors[lower.tri(cors)]))


post.tmb <- as.data.frame(fit)
minsd <- apply(post.tmb, 2, sd) |> min()
maxsd <- apply(post.tmb, 2, sd) |> max()
maxsd/minsd
minsd <- min(sqrt(diag(M)))
maxsd <- max(sqrt(diag(M)))
maxsd/minsd
