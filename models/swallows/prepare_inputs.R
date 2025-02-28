setwd(here::here('models/swallows/'))

data <- readRDS('data.RDS')
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
compile('swallows.cpp')
dyn.load('swallows')
obj <- MakeADFun(data=data, parameters=inits,
                          random=names(inits)[8:10],
                          DLL='swallows')
saveRDS(obj, file='obj.swallows.RDS')


fit <- sample_sparse_tmb(obj, iter=3000, warmup=500,
                         seed=1,
                         control=list(adapt_delta=.95))
pairs_admb(fit, pars=1:5, order='slow')

obj$par |> length()
obj$env$par |> length()
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
